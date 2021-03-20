
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
sys.dt=0.001
sys.temp=3.0
sys.steps=1000
sys.logfreq=10
sys.logfile="test_logger.dat"
sys.logmode="w"
sys.trjdumpfreq=100
#sys.trjdumpfile="test_dump.xyz"
sys.trjdumpfile="test_dump.lammpstrj"
sys.trjdumpmode="w"
sys.box=[9.1456,9.1456,9.1456]
sys.data_file="PS.data"
#sys.thermostat=BerendsenThermostat(0.01)
sys.thermostat=AndersenThermostat(200)
logger=logger_init()
logger.PotentialEnergy=true
logger.KineticEnergy=true
sys.logger=logger
#sys.trajectorydump=XYZTrajectoryDump()
sys.trajectorydump=LAMMPSTrajectoryDump(false,false,true,false)
#sys.bondinteraction=TableBondInteractions(bondtable,bondtable[2,2]-bondtable[1,2])
sys.nonbondinteraction=OriginalhPFInteractions(0.1,[10,10,10])
@time run_hPFMD(sys)