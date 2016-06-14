
--[[
  2D FIVE-MOMENT SIMULATION
  ASYMMETRIC RECONNECTION WITHOUT GUIDE FIELD
  Based on Chen et al. 2016 10.1002/2016GL068243
  NOTES:
  * The L,M,N coordinates in the paper translate to x, y, -z in this script.
--]]

log = function(...) Lucee.logInfo(string.format(...)) end

lightSpeed = 1
mu0 = 1
epsilon0 = 1/math.sqrt(lightSpeed)/mu0

mi = 1
qi = 1
n0 = 1 -- thus wpi0 =1, di0 = 1
mi_me = 25
Ti_Te = 2
n1_n0 = 1/8
wpe_wce = 2
beta0 = 1
beta1 = 1/15
Bz0_B0 = 0
pert0 = 0.36
Bnoise_level = 0.15
Vnoise_level = 0.15
gasGamma = 5/3

me = mi / mi_me
qe = -qi
wpe0 = math.sqrt(mi/me)
wce0 = wpe0/wpe_wce
B0 = -wce0*me/qe
n1 = n0 * n1_n0
B1 = B0 * math.sqrt( (n1/n0) / (beta1/beta0) )
Bz0 = B0 * Bz0_B0
nc = (n0-n1) * ((B0+B1)/2)^2 / (B1^2-B0^2)
pmag0 = B0^2/2/mu0
p0 = pmag0*beta0
kTtotal = p0/n0 -- k_B*(Te+Ti)
p1 = n1*kTtotal
pmag1 = B1^2/2/mu0
vA0 = B0/math.sqrt(mu0*n0*mi)
cs0 = math.sqrt(gasGamma*p0/n0/mi)
vA1 = B1/math.sqrt(mu0*n1*mi)
cs1 = math.sqrt(gasGamma*p1/n1/mi)
wci0 = qi*B0/mi
di0 = 1
de0 = lightSpeed/wpe0
psi0 = pert0 * di0 * B0

l = 1
Lx = 75
Ly = 25
Nx = 1536
Ny = 512
tEnd = 120/wci0
nFrames = 12
Lucee.IsRestarting = false
Lucee.RestartFrame = -1
writeEveryStep = false

cfl = 0.9
limiter = "monotonized-centered"
elcErrorSpeedFactor = 0
mgnErrorSpeedFactor = 1

applyDiff = true
if applyDiff then
   dtRatio = 0.1 -- dtHyp/dtDiff
-- TODO handle nonuniform grid
   dMin = math.min( Lx/Nx, Ly/Ny )
   dtHyp = cfl*dMin/lightSpeed
   alpha = 0.5*dMin^2/dtHyp * dtRatio
   dtDiff = 0.5*dMin^2/alpha
   log("dtHyp/dtDiff = %g", dtHyp/dtDiff)
end

numFluids = 2
charge = {qe, qi}
mass = {me, mi}

writeFwdFields = false
writeBwdFields = false

jzc = 0.5*(B1+B0)/l/mu0
jzc_e = jzc / (1+Ti_Te)
jzc_i = jzc - jzc_e
vzc_e = jzc_e /qe/(nc+0.5*(n1+n0))
vzc_i = jzc_i /qi/(nc+0.5*(n1+n0))

log("====== verifications  ======")
log("B1/B0 = %g", B1/B0)
log("n1/n0 = %g", n1/n0)
log("nc/n0 = %g = 1/%g", nc/n0, n0/nc)
log("vA0/c = %g = 1/%g", vA0/lightSpeed, lightSpeed/vA0)
log("cs0/c = %g = 1/%g", cs0/lightSpeed, lightSpeed/cs0)
log("p0/pmag0 = %g", p0/pmag0)
log("vA1/c = %g = 1/%g", vA1/lightSpeed, lightSpeed/vA1)
log("cs1/c = %g = 1/%g", cs1/lightSpeed, lightSpeed/cs1)
log("p1/pmag1 = %g", p1/pmag1)
log("jzc = %g B0/Ly/mu0", jzc/(B0/Ly/mu0))

log("======= ic values ==========")
log("n0 = %g, n1 = %g, nc+nb = %g", n0, n1, nc+0.5*(n1+n0))
log("p0 = %g, p1 = %g", p0, p1)
log("B0 = %g, B1 = %g", B0, B1)
log("jzc_e = %g", jzc_e)
log("jzc_i = %g", jzc_i)
log("vzc_e = %g = %g vA0", vzc_e, vzc_e/vA0)
log("vzc_i = %g = %g vA0", vzc_i, vzc_i/vA0)

log("======= kinetic parameters =")
log("lightSpeed = %g", lightSpeed)
log("mu0 = %g", mu0)
log("me = %g", me)
log("mi = %g", mi)
log("qe = %g", qe)
log("qi = %g", qi)

log("====== other parameters ====")
log("Bnoise_level = %g", Bnoise_level)
log("Vnoise_level = %g", Vnoise_level)

log("====== normalizations ======")
log("   velocity : %g", vA0)
log("    E field : %g", 2*vA0*B0)
log("temperature : %g", me*vA0^2)
log("       time : %g", 1/wci0)
log("============================")

-----------------------
-- INITIAL CONDITION --
-----------------------
init = function(x,y,z)
   local tanhy = math.tanh(y/l)
   local icosh2y = 1/(math.cosh(y/l))^2
   local Pi = Lucee.Pi

   local nh = nc * icosh2y
   local nb = 0.5*(n1+n0) + 0.5*(n0-n1)*tanhy
   local n = nh+nb
   local rho_e = n*me
   local rho_i = n*mi

   local Bxb = 0.5*(B1-B0) - 0.5*(B1+B0)*tanhy
   local Bx = Bxb + psi0*(Pi/Ly) * math.cos(2*Pi*x/Lx) * math.sin(Pi*y/Ly)
   Bx = Bx * (1 + Bnoise_level*math.random()*math.random(-1,1))
   local By = - psi0 * (2*Pi/Lx) * math.sin(2*Pi*x/Lx) * math.cos(Pi*y/Ly)
   By = By * (1 + Bnoise_level*math.random()*math.random(-1,1))
   local Bz = Bz0

   local rhovx_e, rhovy_e = 0,0
   local rhovx_i, rhovy_i = 0,0
   local jz = 0.5*(B1+B0)/l/mu0 * icosh2y
   local jz_e = jz / (1+Ti_Te) -- current partition due to diamagnetic drift
   local rhovz_e = jz_e * me/qe * (1 + Vnoise_level*math.random()*math.random(-1,1))
   local jz_i = jz - jz_e
   local rhovz_i = jz_i * mi/qi * (1 + Vnoise_level*math.random()*math.random(-1,1))

   local p = n*kTtotal
   local p_e = p / (1+Ti_Te)
   local u_e = p_e / (gasGamma-1) + 0.5*rhovz_e^2/rho_e
   local p_i = p - p_e
   local u_i = p_i / (gasGamma-1) + 0.5*rhovz_i^2/rho_i

   local Ex,Ey,Ez = 0, 0, 0

   return
      rho_e, rhovx_e, rhovy_e, rhovz_e, u_e,
      rho_i, rhovx_i, rhovy_i, rhovz_i, u_i,
      Ex, Ey, Ez, Bx, By, Bz, 0,0
end

----------------------------
-- DECOMPOSITION AND GRID --
----------------------------
decomposition = DecompRegionCalc2D.CartGeneral {}
grid = Grid.RectCart2D {
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {Nx, Ny},
   decomposition = decomposition,
   periodicDirs = {0},
}

------------------------
-- BOUNDARY CONDITION --
------------------------
bcElcCopy = BoundaryCondition.Copy { components = {0, 4} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }

bcIonCopy = BoundaryCondition.Copy { components = {5, 9} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {6, 7, 8} }

bcElcFld = BoundaryCondition.ZeroTangent { components = {10, 11, 12} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {13, 14, 15} }
bcPot = BoundaryCondition.Copy { components = {16, 17}, fact = {-1, 1} }

function createBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      boundaryConditions = {
         bcElcCopy, bcElcWall,
         bcIonCopy, bcIonWall,
         bcElcFld, bcMgnFld, bcPot,
      },
      dir = myDir,
      edge = myEdge,
   }
   return bc
end
bcBottom = createBc(1, "lower")
bcTop = createBc(1, "upper")

function applyBc(myQ, tCurr, tEnd)
   for i,bc in ipairs({bcBottom, bcTop}) do
      bc:setOut( {myQ} )
      bc:advance(tEnd)
   end
   myQ:sync()
end

----------
-- DATA --
----------
createData = function(numComponents)
   if not numComponents then
      numComponents = 18
   end
   return DataStruct.Field2D {
      onGrid = grid,
      numComponents = numComponents,
      ghost = {2, 2},
   }
end
q = createData()
qX = createData()
qNew = createData()
qDup = createData()
if applyDiff then
   qDiff0 = createData(1)
   qDiff = createData(1)
   qJ = createData(1)
   qJ1 = createData(1)
   qJ2 = createData(1)
end

getFields = function(myQ)
   return myQ:alias(0,5), myQ:alias(5,10), myQ:alias(10,18)
end
elc,ion,emf = getFields(q)
elcX,ionX,emfX = getFields(qX)
elcNew,ionNew,emfNew = getFields(qNew)

---------------------------------------
-- HYPERBOLIC EQUATIONS AND SOLVERS --
---------------------------------------
fluidEqn = HyperEquation.Euler { gasGamma = gasGamma, }
fluidEqnLax = HyperEquation.Euler { gasGamma = gasGamma, numericalFlux = "lax" }
emfEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor,
}

createSlvrDir = function(myEqn, input, output, myDir, myLimiter)
   local slvr = Updater.WavePropagation2D {
      onGrid = grid,
      equation = myEqn,
      -- one of no-limiter, zero, min-mod, superbee, 
      -- van-leer, monotonized-centered, beam-warming
      limiter = myLimiter,
      cfl = cfl,
      cflm = 1.1*cfl,
      updateDirections = {myDir}
   }
   slvr:setIn( {input} )
   slvr:setOut( {output} )
   return slvr
end

elcEqnSlvrDir0 = createSlvrDir(fluidEqn, elc, elcX, 0, limiter)
ionEqnSlvrDir0 = createSlvrDir(fluidEqn, ion, ionX, 0, limiter)
emfEqnSlvrDir0 = createSlvrDir(emfEqn, emf, emfX, 0, limiter)

elcEqnSlvrDir1 = createSlvrDir(fluidEqn, elcX, elcNew, 1, limiter)
ionEqnSlvrDir1 = createSlvrDir(fluidEqn, ionX, ionNew, 1, limiter)
emfEqnSlvrDir1 = createSlvrDir(emfEqn, emfX, emfNew, 1, limiter)

elcEqnSlvrDir0Lax = createSlvrDir(fluidEqnLax, elc, elcX, 0, "zero")
ionEqnSlvrDir0Lax = createSlvrDir(fluidEqnLax, ion, ionX, 0, "zero")
emfEqnSlvrDir0Lax = createSlvrDir(emfEqn, emf, emfX, 0, "zero")

elcEqnSlvrDir1Lax = createSlvrDir(fluidEqnLax, elcX, elcNew, 1, "zero")
ionEqnSlvrDir1Lax = createSlvrDir(fluidEqnLax, ionX, ionNew, 1, "zero")
emfEqnSlvrDir1Lax = createSlvrDir(emfEqn, emfX, emfNew, 1, "zero")

----------------------------------
-- HYPERBOLIC EQUATION UPDATERS --
----------------------------------
slvrs = {
   {elcEqnSlvrDir0, ionEqnSlvrDir0, emfEqnSlvrDir0},
   {elcEqnSlvrDir1, ionEqnSlvrDir1, emfEqnSlvrDir1},
}

slvrsLax = {
   {elcEqnSlvrDir0Lax, ionEqnSlvrDir0Lax, emfEqnSlvrDir0Lax},
   {elcEqnSlvrDir1Lax, ionEqnSlvrDir1Lax, emfEqnSlvrDir1Lax},
}

qIn = {q, qX}
qOut = {qX, qNew}

elcOut = {elcX, elcNew}
ionOut = {ionX, ionNew}

function updateHyperEqns(tCurr, tEnd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(tEnd-tCurr)
   local useLaxFlux = false

   for d = 0,1 do
      applyBc(qIn[d+1], tCurr, tEnd)
      for i,slvr in ipairs(slvrs[d+1]) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(tEnd)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end

      if ((fluidEqn:checkInvariantDomain(elcOut[d+1]) == false)
       or (fluidEqn:checkInvariantDomain(ionOut[d+1]) == false)
       or (qOut[d+1]:hasNan())) then
         useLaxFlux = true
      end
   
      if ((myStatus == false) or (useLaxFlux == true)) then
         return myStatus, myDtSuggested, useLaxFlux
      end
   end

   return myStatus, myDtSuggested, useLaxFlux
end

function updateHyperEqnsLax(tCurr, tEnd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(tEnd-tCurr)

   for d = 0,1 do
      applyBc(qIn[d+1], tCurr, tEnd)
      for i,slvr in ipairs(slvrsLax[d+1]) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(tEnd)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end

      if myStatus == false then
         return myStatus, myDtSuggested
      end
   end

   return myStatus, myDtSuggested
end

---------------------
-- SOURCE UPDATERS --
---------------------
srcSlvr = Updater.ImplicitFiveMomentSrc2D {
   onGrid = grid,
   numFluids = numFluids,
   charge = charge,
   mass = mass,
   epsilon0 = epsilon0,
   linearSolver = "partialPivLu",
}

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594–626. doi:10.1016/j.jcp.2013.08.021)
function b(j)
   if (j<2) then 
      return 1.0/3.0
   else 
      return (j^2+j-2)/(2*j*(j+1))
   end
end
--
function a(j) return 1-b(j) end
--
function w1(s) return 4/(s^2+s-2) end
-- 
function mubar(s,j) 
   if (j<2) then 
      return 4/(3*(s^2+s-2)) 
   else 
      return 4*(2*j-1)/(j*(s^2+s-2))*b(j)/b(j-1)
   end
end
--
function mu(j) return (2*j-1)/j*b(j)/b(j-1) end
-- 
function nu(j) return -(j-1)/j*b(j)/b(j-2) end
--
function gbar(s,j) return -a(j-1)*mubar(s,j) end
-- 
function calcNumStages(dhdp) 
   return math.ceil (math.sqrt(4*dhdp+9/4) - 1/2)
end

if applyDiff then
   diffCalc = Updater.RectSecondOrderCentralDiff2D { onGrid = grid }
end

-- function to compute parabolic term
function calcParabolicTerm(qIn, diffOut, t)
   qIn:sync()
   diffCalc:setCurrTime(t)
   diffCalc:setIn( {qIn} )
   diffCalc:setOut( {diffOut} )
   local status, dtSuggested = diffCalc:advance(t)
   diffOut:scale(alpha)
end

function updateSource(elcIn, ionIn, emfIn, tCurr, tEnd)
   local dt = tEnd - tCurr
   srcSlvr:setOut( {elcIn, ionIn, emfIn} )
   srcSlvr:setCurrTime(tCurr)
   srcSlvr:advance(tEnd)

   if applyDiff then
      for dir = 0,2 do
-- applying diffusion on electron only
         local rhov_e = elcIn:alias(dir+1,dir+2)

   local qIn = rhov_e
   local dtRatio = dt/dtDiff
   local numStages = calcNumStages(dtRatio)

   -- we need this in each stage
   calcParabolicTerm(qIn, qDiff0, tCurr)

   -- stage 1
   qJ2:copy(qIn)
   qJ1:combine(1.0, qIn, mubar(numStages,1)*dt, qDiff0)

   -- rest of stages
   for j = 2, numStages do
      calcParabolicTerm(qJ1, qDiff, tCurr)
      qJ:combine(mu(j), qJ1, nu(j), qJ2, 1-mu(j)-nu(j), qIn,
		 mubar(numStages,j)*dt, qDiff, gbar(numStages,j)*dt, qDiff0)

      -- reset fields for next stage
      qJ2:copy(qJ1)
      qJ1:copy(qJ)
   end
   qIn:copy(qJ)

      end
   end
end

----------------------------------------
-- HYPERBOLIC-EQUATION SYSTEM SOLVERS --
----------------------------------------
function updateSystem(tCurr, tEnd)
   local dthalf = 0.5*(tEnd-tCurr)

   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested, useLaxFlux = updateHyperEqns(tCurr, tEnd)

   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)

   return status, dtSuggested, useLaxFlux
end

function updateSystemLax(tCurr, tEnd)
   local dthalf = 0.5*(tEnd-tCurr)

   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested = updateHyperEqnsLax(tCurr, tEnd)

   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)

   return status, dtSuggested
end

------------
-- OUTPUT --
------------
function runOutput(myQ, frame, tCurr, tag, start, stop)
   if not tag then
      tag = "q"
   end
   if start and stop then
      local myQ_ = myQ:alias(start, stop)
      myQ_:write(string.format("%s_%d.h5", tag, frame), tCurr)
   else
      myQ:write(string.format("%s_%d.h5", tag, frame), tCurr)
   end
end

----------
-- MAIN --
----------
function runSimulation(tStart, tEnd, nFrames, initDt)
   local tFrame = (tEnd-tStart)/nFrames
   local tCurr = tStart
   local frame = 1
   local tNextFrame = tCurr + tFrame
   local step = 1
   local stepInFrame = 1
   local myDt = initDt
   local dtSuggested = initDt
   local status = true
   local useLaxFlux = false

   if (Lucee.IsRestarting) then
      rFileName = "q_" .. Lucee.RestartFrame .. ".h5"
      tCurr = q:read(rFileName)
      if not tCurr then
         tCurr = tStart + (tEnd - tStart) * Lucee.RestartFrame / nFrames
      end
      frame = Lucee.RestartFrame + 1
      tNextFrame = tCurr + tFrame
      log('\nRestarting from frame %d tCurr = %g\n', Lucee.RestartFrame, tCurr)
   else
      q:set(init)
   end
   q:sync()
   qNew:copy(q)
   qNew:sync()
   if (not Lucee.IsRestarting) then
      runOutput(qNew, 0, tStart)
   end

   while true do
      if alwaysUseLaxFlux then
         useLaxFlux = true
      end

      qDup:copy(q)

      if (tCurr + myDt > tEnd) then
         myDt = tEnd - tCurr
      end

      if useLaxFlux then
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g; using Lax fluxes",
               step, stepInFrame, frame, tCurr, myDt)
         status, dtSuggested = updateSystemLax(tCurr, tCurr+myDt)
         if status then
            useLaxFlux = false
         end
      else
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g",
               step, stepInFrame, frame, tCurr, myDt)
         status, dtSuggested, useLaxFlux = updateSystem(tCurr, tCurr+myDt)
      end

      if not status then
         log (" ** dt %g too large! Will retake step with dt %g; useLaxFlux = %s",
                            myDt, dtSuggested, tostring(useLaxFlux))
         myDt = dtSuggested
         q:copy(qDup)
      elseif useLaxFlux then
         log(" ** Negative pressure or density! Will retake step with Lax fluxes")
         q:copy(qDup)
      else
         if (qNew:hasNan()) then
            log(" ** NaN occured! Stopping simulation")
            log(" ** Dumping current and last states!")
            runOutput(q, -1, tCurr, "dead")
            runOutput(qNew, 0, tCurr+myDt, "dead")
            break
         end

         q:copy(qNew)
         tCurr = tCurr + myDt

         if (writeEveryStep) then
            log(">>> Writing output at step %d at t = %g...", step, tCurr)
            runOutput(qNew, step, tCurr, "step")
         end

         if writeBwdFields and (tCurr + dtSuggested > tNextFrame) then
            log(">>> Writing output %d (backward) at t = %g...", frame, tCurr)
            runOutput(qNew, frame, tCurr, "b")
         end

         if (tCurr > tNextFrame or tCurr >= tEnd) then
            log(">>> Writing output %d at t = %g...", frame, tCurr)
            runOutput(qNew, frame, tCurr)
            frame = frame + 1
            tNextFrame = tNextFrame + tFrame
            stepInFrame = 0
           
            doWriteFwdFields = true and writeFwdFields
         elseif doWriteFwdFields then
            log(">>> Writing output %d (forward) at t = %g...", frame-1, tCurr)
            runOutput(qNew, frame-1, tCurr, "f")
            doWriteFwdFields = false
         end

         myDt = dtSuggested
         step = step + 1
         stepInFrame = stepInFrame + 1

         if (tCurr >= tEnd) then
            break
         end
      end
   end
end

runSimulation(0, tEnd, nFrames, 100)
