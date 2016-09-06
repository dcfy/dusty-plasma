
--[[
    2D Kink instability
    Harris sheet configuration based on Pritchett 1996
--]]

log = function(...) Lucee.logInfo(string.format(...)) end

lightSpeed = 1
mu0 = 1
epsilon0 = 1/math.sqrt(lightSpeed)/mu0
Pi = math.pi

n0 = 1 -- thus wpi0 = 1, di0 = 1
mi = 1
qi = 1
me = 1/1.0
qe = -1

--physical parameters
--kinetic simulations have an equivalent (w_pe/w_ce = 2.0)
nbOverN0 = 0.2
wpe_over_wce = 2.828282
B0 = 1/wpe_over_wce*math.sqrt(n0*me)

TiOverTe = 1.0
Ttotal = (B0*B0)/2.0/n0
TeFrac = 1.0 / (1.0 + TiOverTe)
TiFrac = 1.0 - TeFrac
Te = Ttotal*TeFrac
Ti = Ttotal*TiFrac
vi = math.sqrt(2*Ti/mi)


vA0 = B0/math.sqrt(n0*mi)
wci0 = B0*qi/mi
larmor_i = vi/wci0

l = 2.0*larmor_i -- Current sheet width
Lx = 25.6*l
Ly = 12.8*l
Nx = 512
Ny = 256
tEnd = 100/wci0
nFrames = 10

cfl = 0.9
limiter = "monotonized-centered"
elcErrorSpeedFactor = 0
mgnErrorSpeedFactor = 1

applyDiff = false
eta_H = 1.5e-5
alpha = eta_H*n0*qe*qe/me
-- to convert etah to a viscosity, multiply by n e^2/m_e
-- TODO make it actually density dependent
dMin = math.min( Lx/Nx, Ly/Ny )
dtHyp = cfl*dMin/lightSpeed
--alpha = 0.5*dMin^2/dtHyp/2.0
dtDiff = 0.5*dMin^2/alpha
canSkipDiff = true
doResetDiff = true

numFluids = 2
charge = {qe, qi}
mass = {me, mi}

writeFwdFields = true
fwdFields_start = 0
fwdFields_stop = 4


di = 1.0
de = di/math.sqrt(mi/me)
elcCollAverageWaveNumber = 1.0/de
ionCollAverageWaveNumber = 1.0/di


--noise generation
Bnoise_level = 0.005
k0 = 1.0 --first wave mode to perturb with noise, 1.0 correspond to box size
kf = 20.0 --last wave mode to perturb with noise
Noise_index = -1.0 --spectral index of the noise
toSeed = true
log("====== verifications  ======")

log("======= ic values ==========")
log("n0 = %g", n0)


log("======= kinetic parameters =")
log("lightSpeed = %g", lightSpeed)
log("mu0 = %g", mu0)
log("me = %g", me)
log("mi = %g", mi)
log("qe = %g", qe)
log("qi = %g", qi)

log("====== other parameters ====")
log("dtHyp/dtDiff = %g", dtHyp/dtDiff)

log("====== normalizations ======")
log("   velocity : %g", vA0)
log("    E field : %g", 2*vA0*B0)
log("temperature : %g", me*vA0^2)
log("       time : %g", 1/wci0)
log("gyroradius_i: %g", vi/wci0)
log("        L   : %g", l)
log("============================")

log("====== other parameters ====")
log("Bnoise_level = %g", Bnoise_level)
log("Noise_index = %g", Noise_index)

-----------------------
--- NOISE FUNCTION ---
-----------------------
noiseGenerator = function(noiseAmp,noiseIndex,kInit,kFinal,x,y) 

   local sin = math.sin
   local cos = math.cos
   local sqrt = math.sqrt
   local BxNoise = 0.0
   local ByNoise = 0.0
   local JzNoise = 0.0

   local xrand_Bx = 0.0
   local xphaserand_Bx = 0.0
   local yrand_Bx = 0.0
   local yphaserand_Bx = 0.0

   local xrand_By = 0.0
   local xphaserand_By = 0.0
   local yrand_By = 0.0
   local yphaserand_By = 0.0
   
   local kindex = (noiseIndex + 1.) / 2.

   if toSeed then
      math.randomseed(x*y*os.time())
      toSeed = false
   end

   for i = math.floor(kInit), math.floor(kFinal) do
      xrand_Bx = math.random(0,2)
      xphaserand_Bx = math.random()
      yrand_Bx = math.random(0,2)
      yphaserand_Bx = math.random()

      xrand_By = math.random(0,2)
      xphaserand_By = math.random()
      yrand_By = math.random(0,2)
      yphaserand_By = math.random()

      BxNoise = BxNoise - xrand_Bx*sin(i*2.0*Pi*y/Ly + 2*Pi*xphaserand_Bx)*i^kindex
      ByNoise = ByNoise + xrand_By*sin(i*2.0*Pi*x/Lx + 2*Pi*xphaserand_By)*i^kindex
      JzNoise = JzNoise + xrand_Bx*(i*2.0*Pi/Ly)*cos(i*2.0*Pi*y/Ly + 2*Pi*xphaserand_Bx)*i^kindex + xrand_By*(i*2.0*Pi/Lx)*cos(i*2.0*Pi*x/Lx + 2*Pi*xphaserand_By)*i^kindex
   end
   BxNoise = noiseAmp*BxNoise
   ByNoise = noiseAmp*ByNoise
   JzNoise = noiseAmp*JzNoise

-- This renormalizes the RMS value to noiseAmp for noiseIndex = -1
   local kdiff = math.floor(kFinal) - math.floor(kInit) + 1.0
   BxNoise = BxNoise/sqrt(2.0*kdiff/3.0)
   ByNoise = ByNoise/sqrt(2.0*kdiff/3.0)
   JzNoise = JzNoise/sqrt(2.0*kdiff/3.0)

   return BxNoise, ByNoise, JzNoise
end


-----------------------
-- INITIAL CONDITION --
-----------------------
init = function(x,y,z)
   local tanh = math.tanh
   local cosh = math.cosh
   local sinh = math.sinh
   local cos = math.cos
   local sin = math.sin
   local sqrt = math.sqrt

   local sech2 = (1.0/cosh(y/l))^2
   local _2pi = 2.0*Pi

   -- background field
   -- compared to the reconnection harris sheet, we send x->-z, z->x
   local Bxb = 0.0
   local Byb = 0.0
   local Bzb = -B0*tanh(y/l)
   -- generate noise
   local BzNoise, ByNoise, JxNoise = noiseGenerator(Bnoise_level, Noise_index, k0, kf, x, y)
   -- addition of perturbation
   local Bx = 0.0
   local By = Byb + ByNoise
   local Bz = Bzb - BzNoise

   local n = n0*nbOverN0 + n0*sech2
   -- pressure balance condition
   local Ttotal = (B0*B0)/2.0/n0
   local TeFrac = 1.0 / (1.0 + TiOverTe)
   local TiFrac = 1.0 - TeFrac
   local Te = Ttotal*TeFrac
   local Ti = Ttotal*TiFrac

   local Jx  = (B0/l)*(-sech2)  + JxNoise
   local Jxe = Jx*TeFrac
   local Jxi = Jx*TiFrac
   local Jze = 0.0
   local Jzi = 0.0

   local rhoe = me*n
   local momze = (me/qe)*Jze
   local momxe = (me/qe)*Jxe
   local uez = momze/rhoe
   local uey = 0.0
   local uex = momxe/rhoe
   local pe = n*Te
   local pezz = n*Te + momze*momze/rhoe
   local pexx = n*Te + momxe*momxe/rhoe
   local pexz = momze*momxe/rhoe
   local peyy = pe
   local peyz = 0.0
   local pexy = 0.0

   local rhoi = mi*n
   local momzi = (mi/qi)*Jzi
   local momxi = (mi/qi)*Jxi
   local pi = n*Ti
   local uiz = momzi/rhoi
   local uiy = 0.0
   local uix = momxi/rhoi
   local pixx = n*Ti + momxi*momxi/rhoi
   local pizz = n*Ti + momzi*momzi/rhoi
   local pixz = momxi*momzi/rhoi
   local pixy = 0.0
   local piyy = pi
   local piyz = 0.0

   return rhoe, momxe, 0.0, momze, pexx, pexy, pexz, peyy, peyz, pezz, 
          rhoi, momxi, 0.0, momzi, pixx, pixy, pixz, piyy, piyz, pizz, 
          0.0, 0.0, 0.0, Bx, By, Bz, 0.0, 0.0
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
-- boundary applicator objects for fluids and fields
bcElcCopy = BoundaryCondition.Copy { components = {0} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcElcPrCopyY = BoundaryCondition.Copy { components = {4, 6, 7, 9} }
bcElcPrFlipY = BoundaryCondition.Copy { components = {5, 8}, fact = {-1, -1} }

bcIonCopy = BoundaryCondition.Copy { components = {10} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {11, 12, 13} }
bcIonPrCopyY = BoundaryCondition.Copy { components = {14, 16, 17, 19} }
bcIonPrFlipY = BoundaryCondition.Copy { components = {15, 18}, fact = {-1, -1} }

bcElcFld = BoundaryCondition.ZeroTangent { components = {20, 21, 22} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {23, 24, 25} }
bcPot = BoundaryCondition.Copy { components = {26, 27} }

-- create boundary condition object
function createBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {
   bcElcCopy, bcElcWall, 
   bcElcPrCopyY, bcElcPrFlipY,
   bcIonCopy, bcIonWall,
   bcIonPrCopyY, bcIonPrFlipY,
   bcElcFld, bcMgnFld, bcPot,
      },
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- create updaters to apply boundary conditions
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
      numComponents = 28
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
   if canSkipDiff then
      rhovDup = createData(3)
   end
end

getFields = function(myQ)
   return myQ:alias(0,10), myQ:alias(10,20), myQ:alias(20,28)
end
elc,ion,emf = getFields(q)
elcX,ionX,emfX = getFields(qX)
elcNew,ionNew,emfNew = getFields(qNew)

---------------------------------------
-- HYPERBOLIC EQUATIONS AND SOLVERS --
---------------------------------------
fluidEqn = HyperEquation.TenMoment { }
fluidEqnLax = HyperEquation.TenMoment { numericalFlux = "lax" }
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
srcSlvr = Updater.ImplicitTenMomentSrc2D {
   onGrid = grid,
   numFluids = numFluids,
   charge = charge,
   mass = mass,
   epsilon0 = epsilon0,
   linearSolver = "analytic",
}

elcCollSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = elcCollAverageWaveNumber,
}
ionCollSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = ionCollAverageWaveNumber,
}

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594626. doi:10.1016/j.jcp.2013.08.021)
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
  
   -- ramp down diffusion to zero at wall boudnaries
   function resetDiff(x,y,z,val)
      return val * 0.5 * (math.cos(2*math.pi*y/Ly) +  1)
   end
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

   -- electron collisional relaxation
   elcCollSrcSlvr:setOut( {elcIn} )
   elcCollSrcSlvr:setCurrTime(tCurr)
   elcCollSrcSlvr:advance(tEnd)

   -- ion collisional relaxation
   ionCollSrcSlvr:setOut( {ionIn} )
   ionCollSrcSlvr:setCurrTime(tCurr)
   ionCollSrcSlvr:advance(tEnd)

   if applyDiff then
      local rhov_e3 = elcIn:alias(1,4)
      if (canSkipDiff) then
         rhovDup:copy(rhov_e3)
      end
      for dir = 0,2 do
-- applying diffusion on electron only
         local rhov_e = elcIn:alias(dir+1,dir+2)
         rhov_e:sync()
         
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

   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested, useLaxFlux = updateHyperEqns(tCurr, tEnd)

   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)

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
