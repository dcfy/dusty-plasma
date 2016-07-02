
--[[
  2D FIVE-MOMENT SIMULATION
  SYMMETRIC RECONNECTION WITHOUT GUIDE FIELD (GEM-LIKE)
  IONS: HYDROGEN+OXYGEN
--]]

log = function(...) Lucee.logInfo(string.format(...)) end

function HERE()
   local info = debug.getinfo(2)
   str = string.format("HERE: %d %s: %s",
       info.currentline, info.source, tostring(info.name))
   log(str)
end

lightSpeed = 1
mu0 = 1
epsilon0 = 1/math.sqrt(lightSpeed)/mu0

n0 = 1
-- H+
-- wpi0 = 1, di0 = 1
mi = 1
qi = 1
-- O+
mo = 16
qo = 1
-- e-
me = 1/25
qe = -1
-- fractions of species temperatures to total temperature
pi_ptot = 1/3
po_ptot = 1/3
pe_ptot = 1/3
-- fractions of species number density to total number density
ni_ntot = 1/2
no_ntot = 1/2
ne_ntot = 1

nb_n0 = 0.2
wpe_wce = 40 -- == vAe/c; needs to be small to make sure vAo/c < 1
beta0 = 1
Bz0_B0 = 0
pert0 = 0.1
Bnoise_level = 0
Vnoise_level = 0
gasGamma = 5/3

nb = n0 * nb_n0
wpe0 = math.sqrt(mi/me)
wce0 = wpe0/wpe_wce
B0 = -wce0*me/qe
Bz0 = B0 * Bz0_B0
pmag0 = B0^2/2/mu0
p0 = pmag0*beta0
Ti0 = p0 * pi_ptot / (n0 * ni_ntot)
To0 = p0 * po_ptot / (n0 * no_ntot)
Te0 = p0 * pe_ptot / (n0 * ne_ntot)
vA0 = B0/math.sqrt(mu0*n0*mi)
cs0 = math.sqrt(gasGamma*p0/n0/mi)
wci0 = qi*B0/mi
di0 = 1
de0 = lightSpeed/wpe0
psi0 = pert0 * di0 * B0

l = 0.9
Lx = 100
Ly = 50
Nx = 512
Ny = 256
tEnd = 200/wci0
nFrames = 20
Lucee.IsRestarting = false
Lucee.RestartFrame = -1

cfl = 0.9
limiter = "monotonized-centered"
elcErrorSpeedFactor = 0
mgnErrorSpeedFactor = 1

applyDiff = true
dtRatio = 0.1 -- dtHyp/dtDiff
-- TODO handle nonuniform grid
dMin = math.min( Lx/Nx, Ly/Ny )
dtHyp = cfl*dMin/lightSpeed
alpha = 0.5*dMin^2/dtHyp * dtRatio
dtDiff = 0.5*dMin^2/alpha
canSkipDiff = true
doResetDiff = true

numFluids = 3
charge = {qe, qi, qo}
mass = {me, mi, mo}

writeFwdFields = false
fwdFields_start = 0
fwdFields_stop = 4

jzc = -B0/l/mu0
jzc_e = jzc * pi_ptot
jzc_i = jzc * pe_ptot
jzc_o = jzc * po_ptot
vzc_e = jzc_e /qe/(n0+nb)
vzc_i = jzc_i /qi/(n0+nb)
vzc_o = jzc_i /qo/(n0+nb)

log("====== verifications  ======")
log("vA0/c = %g = 1/%g", vA0/lightSpeed, lightSpeed/vA0)
log("cs0/c = %g = 1/%g", cs0/lightSpeed, lightSpeed/cs0)
log("p0/pmag0 = %g", p0/pmag0)
log("jzc = %g B0/Ly/mu0", jzc/(B0/Ly/mu0))

log("======= ic values ==========")
log("n0 = %g, nb = %g, n0+nb = %g", n0, nb, n0+nb)
log("p0 = %g", p0)
log("B0 = %g", B0)
log("jzc_e = %g", jzc_e)
log("jzc_i = %g", jzc_i)
log("vzc_e = %g = %g vA0", vzc_e, vzc_e/vA0)
log("vzc_i = %g = %g vA0", vzc_i, vzc_i/vA0)

log("======= kinetic parameters =")
log("lightSpeed = %g", lightSpeed)
log("mu0 = %g", mu0)
log("me = %g", me)
log("mi = %g", mi)
log("mo = %g", mo)
log("qe = %g", qe)
log("qi = %g", qi)
log("qo = %g", qo)

log("====== other parameters ====")
log("dtHyp/dtDiff = %g", dtHyp/dtDiff)

log("====== normalizations ======")
log("============================")

-----------------------
-- INITIAL CONDITION --
-----------------------
init = function(x,y,z)
   local tanhy = math.tanh(y/l)
   local icosh2y = 1/(math.cosh(y/l))^2
   local Pi = Lucee.Pi

   local n = n0 * icosh2y + nb
   local n_e = n * ne_ntot
   local n_i = n * ni_ntot
   local n_o = n * no_ntot
   local rho_e = n_e*me
   local rho_i = n_i*mi
   local rho_o = n_o*mo

   local Bxb = B0*tanhy
   local Bx = Bxb - psi0*(Pi/Ly) * math.cos(2*Pi*x/Lx) * math.sin(Pi*y/Ly) 
   Bx = Bx * (1 + Bnoise_level*math.random()*math.random(-1,1))
   local By = psi0 * (2*Pi/Lx) * math.sin(2*Pi*x/Lx) * math.cos(Pi*y/Ly)
   By = By * (1 + Bnoise_level*math.random()*math.random(-1,1))
   local Bz = Bz0

   local rhovx_e, rhovy_e = 0,0
   local rhovx_i, rhovy_i = 0,0
   local rhovx_o, rhovy_o = 0,0
   local jz = -B0/l/mu0 * icosh2y
   local jz_e = jz * pe_ptot -- current partition due to diamagnetic drift
   local rhovz_e = jz_e * me/qe * (1 + Vnoise_level*math.random()*math.random(-1,1))
   local jz_i = jz * pi_ptot
   local rhovz_i = jz_i * mi/qi * (1 + Vnoise_level*math.random()*math.random(-1,1))
   local jz_o = jz * po_ptot
   local rhovz_o = jz_o * mo/qo * (1 + Vnoise_level*math.random()*math.random(-1,1))

   local p_e = n * Te0
   local u_e = p_e / (gasGamma-1) + 0.5*rhovz_e^2/rho_e
   local p_i = n * Ti0
   local u_i = p_i / (gasGamma-1) + 0.5*rhovz_i^2/rho_i
   local p_o = n * To0 
   local u_o = p_o / (gasGamma-1) + 0.5*rhovz_o^2/rho_o

   local Ex,Ey,Ez = 0, 0, 0

   return
      rho_e, rhovx_e, rhovy_e, rhovz_e, u_e,
      rho_i, rhovx_i, rhovy_i, rhovz_i, u_i,
      rho_o, rhovx_o, rhovy_o, rhovz_o, u_o,
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

bcOxyCopy = BoundaryCondition.Copy { components = {10, 14} }
bcOxyWall = BoundaryCondition.ZeroNormal { components = {11, 12, 13} }

bcElcFld = BoundaryCondition.ZeroTangent { components = {15, 16, 17} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {18, 19, 20} }
bcPot = BoundaryCondition.Copy { components = {21, 22}, fact = {-1, 1} }

function createBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      boundaryConditions = {
         bcElcCopy, bcElcWall,
         bcIonCopy, bcIonWall,
         bcOxyCopy, bcOxyWall,
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
      numComponents = numFluids*5+8
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
   qDiff = createData(1)
   if canSkipDiff then
      rhovDup = createData(3)
   end
end

getFields = function(myQ)
   return myQ:alias(0,5), myQ:alias(5,10), myQ:alias(10,15), myQ:alias(15,23)
end
elc,ion,oxy,emf = getFields(q)
elcX,ionX,oxyX,emfX = getFields(qX)
elcNew,ionNew,oxyNew,emfNew = getFields(qNew)

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
oxyEqnSlvrDir0 = createSlvrDir(fluidEqn, oxy, oxyX, 0, limiter)
emfEqnSlvrDir0 = createSlvrDir(emfEqn, emf, emfX, 0, limiter)

elcEqnSlvrDir1 = createSlvrDir(fluidEqn, elcX, elcNew, 1, limiter)
ionEqnSlvrDir1 = createSlvrDir(fluidEqn, ionX, ionNew, 1, limiter)
oxyEqnSlvrDir1 = createSlvrDir(fluidEqn, oxyX, oxyNew, 1, limiter)
emfEqnSlvrDir1 = createSlvrDir(emfEqn, emfX, emfNew, 1, limiter)

elcEqnSlvrDir0Lax = createSlvrDir(fluidEqnLax, elc, elcX, 0, "zero")
ionEqnSlvrDir0Lax = createSlvrDir(fluidEqnLax, ion, ionX, 0, "zero")
oxyEqnSlvrDir0Lax = createSlvrDir(fluidEqnLax, oxy, oxyX, 0, "zero")
emfEqnSlvrDir0Lax = createSlvrDir(emfEqn, emf, emfX, 0, "zero")

elcEqnSlvrDir1Lax = createSlvrDir(fluidEqnLax, elcX, elcNew, 1, "zero")
ionEqnSlvrDir1Lax = createSlvrDir(fluidEqnLax, ionX, ionNew, 1, "zero")
oxyEqnSlvrDir1Lax = createSlvrDir(fluidEqnLax, oxyX, oxyNew, 1, "zero")
emfEqnSlvrDir1Lax = createSlvrDir(emfEqn, emfX, emfNew, 1, "zero")

----------------------------------
-- HYPERBOLIC EQUATION UPDATERS --
----------------------------------
slvrs = {
   {elcEqnSlvrDir0, ionEqnSlvrDir0, oxyEqnSlvrDir0, emfEqnSlvrDir0},
   {elcEqnSlvrDir1, ionEqnSlvrDir1, oxyEqnSlvrDir1, emfEqnSlvrDir1},
}

slvrsLax = {
   {elcEqnSlvrDir0Lax, ionEqnSlvrDir0Lax, oxyEqnSlvrDir0Lax, emfEqnSlvrDir0Lax},
   {elcEqnSlvrDir1Lax, ionEqnSlvrDir1Lax, oxyEqnSlvrDir1Lax, emfEqnSlvrDir1Lax},
}

qIn = {q, qX}
qOut = {qX, qNew}

elcOut = {elcX, elcNew}
ionOut = {ionX, ionNew}
oxyOut = {oxyX, oxyNew}

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
       or (fluidEqn:checkInvariantDomain(oxyOut[d+1]) == false)
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

if applyDiff then
   diffCalc = Updater.RectSecondOrderCentralDiff2D { onGrid = grid }
  
   -- ramp down diffusion to zero at wall boudnaries
   function resetDiff(x,y,z,val)
      return val * 0.5 * (math.cos(2*math.pi*y/Ly) +  1)
   end

end

function updateSource(elcIn, ionIn, oxyIn, emfIn, tCurr, tEnd)
   local dt = tEnd - tCurr
   srcSlvr:setOut( {elcIn, ionIn, oxyIn, emfIn} )
   srcSlvr:setCurrTime(tCurr)
   srcSlvr:advance(tEnd)

   if applyDiff then
      local rhov_e3 = elcIn:alias(1,4)
      if (canSkipDiff) then
         rhovDup:copy(rhov_e3)
      end
      for dir = 0,2 do
-- applying diffusion on electron only
         local rhov_e = elcIn:alias(dir+1,dir+2)
         rhov_e:sync()

         diffCalc:setCurrTime(tCurr)
         diffCalc:setIn( {rhov_e} )
         diffCalc:setOut( {qDiff} )
         diffCalc:advance(tEnd)
         if doResetDiff then
            qDiff:set(resetDiff)
         end

         rhov_e:accumulate(alpha*dt, qDiff)
      end
      if (canSkipDiff and fluidEqn:checkInvariantDomain(elcIn) == false) then
         log(" ** Parabolic source leads to negative pressure. Will skip it.")
         rhov_e3:copy(rhovDup)
      end
   end
end

----------------------------------------
-- HYPERBOLIC-EQUATION SYSTEM SOLVERS --
----------------------------------------
function updateSystem(tCurr, tEnd)
   local dthalf = 0.5*(tEnd-tCurr)

   -- check if time-step is too large for stability
   if applyDiff then
      local cflm = 0.51
      local cfl = 0.5
      local cfla = alpha*dthalf/dMin^2
      if (cfla > cflm) then
         return false, dthalf*cfl/cfla
      end
   end

   updateSource(elc, ion, oxy, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested, useLaxFlux = updateHyperEqns(tCurr, tEnd)

   updateSource(elcNew, ionNew, oxyNew, emfNew, tCurr, tCurr+dthalf)

   return status, dtSuggested, useLaxFlux
end

function updateSystemLax(tCurr, tEnd)
   local dthalf = 0.5*(tEnd-tCurr)

   -- check if time-step is too large for stability
   if applyDiff then
      local cflm = 0.51
      local cfl = 0.5
      local cfla = alpha*dthalf/dMin^2
      if (cfla > cflm) then
         return false, dthalf*cfl/cfla
      end
   end

   updateSource(elc, ion, oxy, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested = updateHyperEqnsLax(tCurr, tEnd)

   updateSource(elcNew, ionNew, oxyNew, emfNew, tCurr, tCurr+dthalf)

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
            break
         end

         q:copy(qNew)
         tCurr = tCurr + myDt
         if (tCurr > tNextFrame or tCurr >= tEnd or outputEveryStep) then
            log(">>> Writing output %d at t = %g...", frame, tCurr)
            runOutput(qNew, frame, tCurr)
            frame = frame + 1
            tNextFrame = tNextFrame + tFrame
            stepInFrame = 0
           
            doWriteFwdFields = true and writeFwdFields
         end
         if doWriteFieldsFwd then
            log(">>> Writing output %d (forward) at t = %g...", frame-1, tCurr)
            runOutput(qNew, frame-1, tCurr, "f", fwdFields_start, fwdFields_stop)
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
