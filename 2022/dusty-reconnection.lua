-- 5-moment dust plasma modeling --
---------------------------
-- PREAMBLE; DONT CHANGE --
---------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local Logger = require "Logger"
​
local logger = Logger {logToFile = True}
​
local log = function(...)
   logger(string.format(...))
   logger("\n")
end
​
local sqrt, sin, cos = math.sqrt, math.sin, math.cos
local tanh, cosh = math.tanh, math.cosh
​
---------------
-- PARMETERS --
---------------
local gamma = 5.0 / 3.0 -- gas gamma
local mu0 = 1.0
local c = 1.0 -- light speed
local q = 1.0 -- unit charge
local me = 1.0 -- electron mass
local ne0 = 1.0 -- reference electron number density
​
local beta = 1.0 -- plasma beta, p / pmag
local nb__n0 = 0.2 -- background density
local delta__dd0 = 0.5 -- used to determine transition layer thickness delta
local pert = 1e-3 -- relative perturbation level, psi0 / (B0 * dd0)
​
local vAe0__c = 0.5 -- used to determine vAe0, then b0
local mi__me = 25.0 -- mi/me
local Zd = 10.0 -- dust charge number, qd = -Zd * q
local di0__de0 = 4.0 -- dust-electron inertial length ratio
local dd0__di0 = 5.0 -- dust-ion inertial length ratio
local Ti0__Te0 = 1.0 -- Ti/Te
local Td0__Te0 = 1.0 -- Td/Te
​
local eps0 = 1 / c ^ 2 / mu0
local qi = q
local qe = -q
local qd = -Zd * q
local mi = me * mi__me
local de0 = sqrt(me / (ne0 * qe ^ 2 * mu0)) -- refernce electron inertial length
local di0 = de0 * di0__de0 -- refernce ion inertial length
local ni0 = mi / qi ^ 2 / mu0 / di0 ^ 2 -- reference ion number density
local nd0 = (ni0 - ne0) / Zd -- reference dusty number density
local dd0 = di0 * dd0__di0 -- reference dusty inertial length
local md = nd0 * qd ^ 2 * mu0 * dd0 ^ 2 -- dusty mass
local vAe0 = c * vAe0__c -- reference electron Alfven speed
local B0 = vAe0 * sqrt(ne0 * me) -- reference B field strength
local pmag0 = B0 ^ 2 / (2 * mu0) -- reference magnetic pressure
local p0 = pmag0 * beta -- reference pressure
-- uniform electron, ion, and dust temperature
local Te0 = p0 / (ne0 + nd0 * Td0__Te0 + ni0 * Ti0__Te0)
local Ti0 = Te0 * Ti0__Te0
local Td0 = Te0 * Td0__Te0
local delta = dd0 * delta__dd0 -- transition layer thickness
local psi0 = B0 * dd0 * pert -- perturbation flux function magnitude
local wpe0, wpi0, wpd0 = c / de0, c / di0, c / dd0 -- plasma frequencies
-- reference cyclotron frequencies
local wci0 = math.abs(B0 / mi / qi)
local wce0 = math.abs(B0 / me / qe)
local wcd0 = math.abs(B0 / md / qd)
-- reference thermal speed
local vti0 = sqrt(Ti0 / mi)
local vte0 = sqrt(Te0 / me)
local vtd0 = sqrt(Td0 / md)
​
local Lx, Ly = 51.2 * dd0, 25.6 * dd0 -- domain size
local Nx, Ny = 1600, 800 -- grid  size
-- # of cpus in each direction; or nil for automatic assignment
local decompCuts = {40, 40}
local cfl = 0.75 -- CFL number; set to smaller than 0.95
local tEnd = 25.0 / wcd0 -- end simulation at this time
local nFrame = 250 -- number of output frames
local tFrame = tEnd / nFrame
local frames_between_restart = 5 -- write a restart file every xxx frames
local restartFrameEvery = frames_between_restart / nFrame
local dx, dy = Lx / Nx, Ly / Ny
local dMin = math.min(dx, dy)
local dt = dMin / c * cfl
​
log("%45s = %g", "gamma", gamma)
log("%45s = %g", "mu0", mu0)
log("%45s = %g", "lightSpeed", c)
log("%45s = %g", "epsilon0", eps0)
log("%45s = %g", "plasma beta", beta)
log("%45s = %g", "transition layer thickness: delta", delta)
log("%45s = %g, %g, %g", "mi, me, md", mi, me, md)
log("%45s = %g, %g, %g", "mi/me, md/me, md/mi", mi / me, md / me, md / mi)
log("%45s = %g", "Zd", Zd)
log("%45s = %g, %g, %g", "qi, qe, qd", qi, qe, qd)
log("%45s = %g, %g, %g", "ni0, ne0, nd0", ni0, ne0, nd0)
log("%45s = %g, %g, %g", "inertial lengths: di0, de0, dd0", di0, de0, dd0)
log("%45s = %g, %g", "initial temperature ratios: Ti0/Te0, Td0/Te0",
   Ti0__Te0, Td0__Te0)
log("%45s = %g, %g, %g", "vti0, vte0, vtd0", vti0, vte0, vtd0)
log("%45s = %g, %g, %g", "wpi0, wpe0, wpd0", wpi0, wpe0, wpd0)
log("%45s = %g, %g, %g", "wci0, wce0, wcd0", wci0, wce0, wcd0)
log("%45s = %g = %g di0 = %g de0 = %g dd0", "Lx", Lx, Lx / di0, Lx / de0,
    Lx / dd0)
log("%45s = %g = %g di0 = %g de0 = %g dd0", "Ly", Ly, Ly / di0, Ly / de0,
    Ly / dd0)
log("%45s = %g, %g", "Nx, Ny", Nx, Ny)
log("%45s = %g, %g", "dx, dy", dx, dy)
log("%45s = %g = delta / %g", "dMin = min(dx, dy)", dMin, delta / dMin)
log("%45s = di0 / %g = de0 / %g = dd0 / %g", "", di0 / dMin, de0 / dMin, dd0 / dMin)
log("%45s = %g = %g/wci0 = %g/wcd0", "tEnd", tEnd, tEnd * wci0, tEnd * wcd0)
log("%45s = %g", "nFrame", nFrame)
log("%45s = %g = %g/wci0 = %g/wcd0", "tEnd/nFrame", tFrame, tFrame * wci0,
   tFrame * wcd0)
log("%45s = %g", "restartFrameEvery", restartFrameEvery)
log("%45s = %g", "cfl", cfl)
log("%45s = %g = %g / wpe0 = %g / wpi0", "dt = dMin / lightSpeed", dt,
    dt * wpe0, dt * wpi0)
log("%45s = %g", "estimated # of steps", math.ceil(tEnd / dt))
​
---------------------------------
-- CREATING THE SIMULATION APP --
---------------------------------
local momentApp = Moments.App {
   logToFile = true,
​
   timeStepper = "fvDimSplit",
   lower = {-Lx / 2, -Ly / 2},
   upper = {Lx / 2, Ly / 2},
   cells = {Nx, Ny},
   decompCuts = decompCuts,
   periodicDirs = {1}, -- periodic directions, 1: x, 2: y, 3: z
   cflFrac = cfl,
   suggestedDt = dt,
   tEnd = tEnd,
   nFrame = nFrame,
   restartFrameEvery = restartFrameEvery,
​
   -- electron species
   elc = Moments.Species {
      charge = qe,
      mass = me,
      equation = Euler {gasGamma = gamma},
      equationInv = Euler {gasGamma = gamma, numericalFlux = "lax"},
      init = function(t, xn)
         local x, y = xn[1], xn[2]
         local sech2 = (1.0 / math.cosh(y / delta)) ^ 2
         local ne = ne0 * (sech2 + nb__n0)
         local rhoe = ne * me
         local rhovxe = 0.0
         local rhovye = 0.0
         local Jz = -(B0 / delta) * sech2
         local Jze = Jz * Te0 / (Te0 + Td0 + Ti0)
         local rhovze = Jze * me / qe
         local energye = ne * Te0 / (gamma - 1) + 0.5 * rhovze * rhovze / rhoe
         return rhoe, rhovxe, rhovye, rhovze, energye
      end,
      bcy = {Euler.bcWall, Euler.bcWall}
   },
​
   -- dust species
   dust = Moments.Species {
      charge = qd,
      mass = md,
      equation = Euler {gasGamma = gamma},
      equationInv = Euler {gasGamma = gamma, numericalFlux = "lax"},
      init = function(t, xn)
         local x, y = xn[1], xn[2]
         local sech2 = (1.0 / math.cosh(y / delta)) ^ 2
         local nd = nd0 * (sech2 + nb__n0)
         local rhod = nd * md
         local rhovxd = 0.0
         local rhovyd = 0.0
         local Jz = -(B0 / delta) * sech2
         local Jzd = Jz * Ti0 / (Te0 + Td0 + Ti0)
         local rhovzd = Jzd * md / qd
         local energyd = nd * Te0 / (gamma - 1) + 0.5 * rhovzd * rhovzd / rhod
         return rhod, rhovxd, rhovyd, rhovzd, energyd
      end,
      bcy = {Euler.bcWall, Euler.bcWall}
   },
​
   -- ion species
   ion = Moments.Species {
      charge = qi,
      mass = mi,
      equation = Euler {gasGamma = gamma},
      equationInv = Euler {gasGamma = gamma, numericalFlux = "lax"},
      init = function(t, xn)
         local x, y = xn[1], xn[2]
         local sech2 = (1.0 / math.cosh(y / delta)) ^ 2
         local ni = ni0 * (sech2 + nb__n0)
         local rhoi = ni * mi
         local rhovxi = 0.0
         local rhovyi = 0.0
         local Jz = -(B0 / delta) * sech2
         local Jzi = Jz * Ti0 / (Te0 + Td0 + Ti0)
         local rhovzi = Jzi * mi / qi
         local energyi = ni * Te0 / (gamma - 1) + 0.5 * rhovzi * rhovzi / rhoi
         return rhoi, rhovxi, rhovyi, rhovzi, energyi
      end,
      bcy = {Euler.bcWall, Euler.bcWall}
   },
​
   field = Moments.Field {
      epsilon0 = eps0,
      mu0 = mu0,
      init = function(t, xn)
         local x, y = xn[1], xn[2]
         local PI = math.pi
         local Bx0 = B0 * tanh(y / delta)
         local By0 = 0.0
         local Bx1 = -psi0 * (PI / Ly) * cos(2 * PI * x / Lx) * sin(PI * y / Ly)
         local By1 = psi0 * (2 * PI / Lx) * sin(2 * PI * x / Lx) *
                         cos(PI * y / Ly)
         local Bx = Bx0 - Bx1
         local By = By0 - By1
         local Bz = 0.0
         return 0.0, 0.0, 0.0, Bx, By, Bz
      end,
      bcy = {Moments.Field.bcReflect, Moments.Field.bcReflect}
   },
​
   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "dust", "ion"},
      timeStepper = "direct"
   }
​
}
​
--------------------------------
-- RUNNING THE SIMULATION APP --
--------------------------------
momentApp:run()
