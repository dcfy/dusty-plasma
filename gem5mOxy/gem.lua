
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
epsilon0 = 1/lightSpeed^2/mu0
gasGamma = 5/3

-- mass and charge of each species
-- e-
me = 1/25
qe = -1
-- H+
mi = 1
qi = 1
-- O+
mo = 16
qo = 1
-- fractions of species pressures to total pressure
-- used to decompose current according to diamagnetic drift
pe_frac = 1/3
pi_frac = 1/3
po_frac = 1/3
-- fractions of species number density to total number density
-- negative charged species fractions should add up to 1
-- positive charged species fractions should add up to 1
ne_frac = 1
ni_frac = 1/2
no_frac = 1/2

-- characteristic number density of all negatively/positively charged species
n0 = 1
-- to compute background number density nb
nb_n0 = 0.2
-- == vAe/c, used to determine B0
wpe0_wce0 = 10
-- p0/B0, p0 being the characteristic total pressure
beta0 = 1
-- guide field
Bz0_B0 = 0
-- initial perturbation to form island; see psi0
pert0 = 0.1
-- random noise to break symmetry if necessary
Bnoise_level = 0
Vnoise_level = 0

-- derived parameters
-- background number density
nb = n0 * nb_n0
-- characteristic B field
B0 = lightSpeed * math.sqrt(mu0 * n0 * ne_frac * me) / wpe0_wce0
-- guide field
Bz0 = B0 * Bz0_B0
-- magnetic pressured based on B0
pmag0 = B0^2/2/mu0
-- characteristic total pressure
p0 = pmag0*beta0
-- uniform temperatures of each species
Te0 = p0 * pe_frac / (n0 * ne_frac)
Ti0 = p0 * pi_frac / (n0 * ni_frac)
To0 = p0 * po_frac / (n0 * no_frac)

wci0 = qi*B0/mi
-- hydrogen plasma frequency based on number density n0*ni_frac
wpi0 = math.sqrt(n0 * ni_frac * qi^2 / epsilon0 / mi)
-- hydrogen inertia length based on wpi0
di0 = lightSpeed/wpi0

-- perturbation magnetic flux function
psi0 = pert0 * di0 * B0
-- initial current sheet thickness
l = 0.9*di0

-- box size
Lx = 100*di0
Ly = 50*di0
-- grid size
Nx = 512
Ny = 256

-- run control
tEnd = 200/wci0
nFrames = 20
Lucee.IsRestarting = false
Lucee.RestartFrame = -1

-- other numerical parameters
cfl = 0.9
limiter = "monotonized-centered"
elcErrorSpeedFactor = 0
mgnErrorSpeedFactor = 1

-- additional diffusion (hyper-resistivity) term to electron momentum equation
applyDiff = true
-- dtHyp/dtDiff, larger value means larger diffusion and smaller time step size
dtRatio = 0.1
-- smallest grid size, used to compute time step size
dMin = math.min( Lx/Nx, Ly/Ny )
-- time step size due to cfl condition of hyperbolic equation
dtHyp = cfl*dMin/lightSpeed
alpha = 0.5*dMin^2/dtHyp * dtRatio -- diffusion coefficient
dtDiff = 0.5*dMin^2/alpha -- should match dtHyp / dtRatio
canSkipDiff = true -- skip diffusion term if negative pressure/density occurs
doResetDiff = true -- reset diffusion term, e.g., set it to zero near non-periodic boundaries

numFluids = 3
charge = {qe, qi, qo}
mass = {me, mi, mo}

-- switch to write all data in the step next to a main i/o frame
-- e.g., if nFrame = 10, tEnd = 1, then the first main frame writes at tEnd = 0.1
-- assume the step count is, say, 1000, then a fwd (forward) field will be written
-- at step 1001
writeFwdFields = false
-- first component of fields to be written
fwdFields_start = 0 -- electron density
-- last component of fields to be written + 1
fwdFields_stop = 5 -- electron energy + 1
-- .e.g., fwdFields_start, stop = 0, 5 will write all electron fields

-- diagnostic parameters
vAe0 = B0/math.sqrt(mu0*n0*me)
vAi0 = B0/math.sqrt(mu0*n0*mi)
vAo0 = B0/math.sqrt(mu0*n0*mo)
cse0 = math.sqrt(gasGamma*p0/n0/me)
csi0 = math.sqrt(gasGamma*p0/n0/mi)
cso0 = math.sqrt(gasGamma*p0/n0/mo)

-- J=curl B
jzc = -B0/l/mu0
-- decomposition of current due to diamagnetic drift
jzc_e = jzc * pi_frac
jzc_i = jzc * pe_frac
jzc_o = jzc * po_frac
-- drift velocities
vzc_e = jzc_e /qe/(n0+nb)
vzc_i = jzc_i /qi/(n0+nb)
vzc_o = jzc_i /qo/(n0+nb)

log("====== verifications  ======")
log("vAe0/c = %g = 1/%g", vAe0/lightSpeed, lightSpeed/vAe0)
log("vAi0/c = %g = 1/%g", vAi0/lightSpeed, lightSpeed/vAi0)
log("vAo0/c = %g = 1/%g", vAo0/lightSpeed, lightSpeed/vAo0)
log("csi0/c = %g = 1/%g", csi0/lightSpeed, lightSpeed/csi0)
log("cse0/c = %g = 1/%g", cse0/lightSpeed, lightSpeed/cse0)
log("cso0/c = %g = 1/%g", cso0/lightSpeed, lightSpeed/cso0)
log("p0/pmag0 = %g", p0/pmag0)
log("jzc = %g B0/Ly/mu0", jzc/(B0/Ly/mu0))

log("======= ic values ==========")
log("n0 = %g, nb = %g, n0+nb = %g", n0, nb, n0+nb)
log("p0 = %g", p0)
log("B0 = %g", B0)
log("jzc_e = %g", jzc_e)
log("jzc_i = %g", jzc_i)
log("jzc_o = %g", jzc_o)
log("vzc_e = %g = %g vAi0", vzc_e, vzc_e/vAi0)
log("vzc_i = %g = %g vAi0", vzc_i, vzc_i/vAi0)
log("vzc_o = %g = %g vAi0", vzc_o, vzc_o/vAi0)

log("======= kinetic parameters =")
log("lightSpeed = %g", lightSpeed)
log("mu0 = %g", mu0)
log(" ====== masses")
log("me = %g", me)
log("mi = %g", mi)
log("mo = %g", mo)
log(" ====== charges")
log("qe = %g", qe)
log("qi = %g", qi)
log("qo = %g", qo)
log(" ====== pressure fractions")
log("pe_frac = %g", pe_frac)
log("pi_frac = %g", pi_frac)
log("po_frac = %g", po_frac)
log(" ====== number density fractions")
log("ne_frac = %g", ne_frac)
log("ni_frac = %g", ni_frac)
log("no_frac = %g", no_frac)

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
   local n_e = n * ne_frac
   local n_i = n * ni_frac
   local n_o = n * no_frac
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
   local jz_e = jz * pe_frac -- current partition due to diamagnetic drift
   local rhovz_e = jz_e * me/qe * (1 + Vnoise_level*math.random()*math.random(-1,1))
   local jz_i = jz * pi_frac
   local rhovz_i = jz_i * mi/qi * (1 + Vnoise_level*math.random()*math.random(-1,1))
   local jz_o = jz * po_frac
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
-- create arrays to store data
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
-- input of sweep along x
q = createData()
-- output of sweep along x and input of sweep along y
qX = createData()
-- output of sweep along y and input of sweep along z
qNew = createData()
-- used to save old data in case the whole step should be re-taken
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
-- default equation using Roe flux
fluidEqn = HyperEquation.Euler { gasGamma = gasGamma, }
-- posistivity-conserving equation using Lax flux
fluidEqnLax = HyperEquation.Euler { gasGamma = gasGamma, numericalFlux = "lax" }
-- perfectly hyperbolic Maxwell equation
emfEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   -- divE and divB cleaning
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor,
}

-- function to get solver of given equation (myEqn) along a given direction (myDir)
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

-- myStatus can be false if time step size is too large
-- useLaxFlux can be true if negative density/pressure occurs
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

-- apply diffusion on momentum equation(s)
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
         -- duplicate all momentum terms, would be copied back if 
         -- diffusion causes negative density/pressure
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
         -- status can be false if time step size is too large,
         -- will retake with new, smaller time step size
      else
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g",
               step, stepInFrame, frame, tCurr, myDt)
         status, dtSuggested, useLaxFlux = updateSystem(tCurr, tCurr+myDt)
         -- status is false if time step size is too large,
         -- will retake with new, smaller time step size;
         -- useLaxFlux is true if negative density/pressure occurs,
         -- will retake with Lax flux
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
