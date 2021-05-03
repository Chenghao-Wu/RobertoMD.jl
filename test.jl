
include("./hPF-MD/hPF-MD.jl")

using DelimitedFiles
using .hPFMD

# read bondtable forcefield
function read_bondtable(filename::AbstractString)
    bondtable=readdlm(filename, '\t', Float32, '\n')
    return bondtable
end

bondtable=read_bondtable("bondAA.pot.tablenew")

sys=sys_init()
sys.initial_velocity="zero"
sys.dt=0.06
sys.temp=1.0
sys.steps=100000
sys.logfreq=100
sys.logfile="test_logger_hPFMD1.dat"
sys.logmode="w"
sys.trjdumpfreq=10000
#sys.trjdumpfile="test_dump.xyz"
sys.trjdumpfile="test_dump_hPFMD1.lammpstrj"
sys.trjdumpmode="w"
sys.box=[11.0,11.0,11.0]
sys.data_file="input2.data"
sys.thermostat=BerendsenThermostat(1)
#sys.thermostat=AndersenThermostat(2)
logger=logger_init()
logger.PotentialEnergy=true
logger.KineticEnergy=true
logger.Momentum=true
sys.logger=logger
#sys.trajectorydump=XYZTrajectoryDump()
sys.trajectorydump=LAMMPSTrajectoryDump(false,false,true,false)
#sys.bondinteraction=TableBondInteractions(bondtable,bondtable[2,2]-bondtable[1,2])
sys.bondinteraction=HarmonicBondInteractions(100,0.5)
mesh=init_mesh([11.0,11.0,11.0],10,10,10)
sys.nonbondinteraction=originalhPFInteractions(mesh,0.2)
#sys.nonbondinteraction=spectralhPFInteractions(mesh,0.2,0.5) # Mesh, kappa, sigma
@time run_hPFMD(sys)