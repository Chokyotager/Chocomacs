import os
import subprocess
import shutil
import json
import math

from logger import Logger
import mdpreader

config = json.loads(open("config.json").read())

def execute_gmx (command, step=1, precommand=None):

    assert type(step) == int

    command = command.replace("$OUT_DIR$", config["output-dir"] + "/")
    command = command.replace("$STEP$", str(step))
    command = command.replace("$LAST_STEP$", str(step - 1))
    command = config["gromacs-bin"] + " " + command

    if precommand != None:

        command = precommand + " | " + command

    output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if output.returncode != 0:
        logger.log(4, output.stdout)
        logger.log(4, "Crashed at command: {}".format(command))
        logger.log(4, "Chocomacs has crashed. Forcibly exiting.")
        exit()

    logger.log(1, output.stdout)

    return output

def get_parameter (output, characteristic):

    stdout = output.stdout.split("\n")

    for output_line in stdout:

        if output_line.startswith(characteristic):

            return output_line.replace(characteristic, str()).split()

os.system("rm -rf test_out")

output_directory = config["output-dir"]

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

if not os.path.exists(output_directory + "/input"):
    os.makedirs(output_directory + "/input")

logger = Logger(log_level=config["log-level"], log_file=output_directory + "/log.txt")

logger.log(2, "\n\n\n█▀▀ █▄█ █▀█ █▀▀ █▀█ █▄ ▄█ ▄▀▄ █▀▀ █▀▀ \n█▄▄ █ █ █▄█ █▄▄ █▄█ █ ▀ █ █▀█ █▄▄ ▄██ ")
logger.log(2, "\n[Version 1.0, by Hilbert Lam]\n\n")

# Start by generating the FF / combining topology file
# Self-make an index from the GRO file

logger.log(2, "Running run name {}".format(config["run-name"]))

run_file_name = config["run-file"].split("/")[-1]

topology = config["topology"]

input_file = output_directory + "/input/" + run_file_name

current_step = 1

# Generate FF
shutil.copy("config.json", output_directory + "/input/config.json")
shutil.copy(config["run-file"], input_file)

logger.log(2, "Generating topology")

execute_gmx("pdb2gmx -f {} -o $OUT_DIR$/step$STEP$_read.gro -i $OUT_DIR$/posre.itp -p $OUT_DIR$/topol.top -n $OUT_DIR$/index.ndx -ff {} -water {} -ignh -ter".format(input_file, topology["force-field"], topology["water-model"]), step=current_step)

current_step += 1

logger.log(2, "Generating periodic boundaries")

# Centre and PBC accordingly
precommand = "printf '%s\n' 0"
output = execute_gmx("editconf -f $OUT_DIR$/step$LAST_STEP$_read.gro -princ -o $OUT_DIR$/step$STEP$_pbc.gro -bt {} -d {}".format(config["pbc"], str(config["pbc-boundary-nm"])), step=current_step, precommand=precommand)
last_out = config["output-dir"] + "/step{}_pbc.gro".format(current_step)

current_step += 1

box_size = [float(x) for x in get_parameter(output, "new system size :")]

logger.log(2, "Box size is: {} (nm)".format(box_size))

solvation = config["solvation"]

if solvation["water"]:

    logger.log(2, "Performing solvation")

    # Solvate
    execute_gmx("solvate -cp $OUT_DIR$/step$LAST_STEP$_pbc.gro -o $OUT_DIR$/step$STEP$_solvation.gro -p $OUT_DIR$/topol.top", step=current_step)

    last_out = "step{}_solvation.gro".format(current_step)

    current_step += 1

    ions = solvation["ions"]

    if len(ions) > 0:

        logger.log(2, "Adding ions")

        shutil.copy("templates/arbitrary.mdp", config["output-dir"] + "/arbitrary.mdp")
        shutil.copy(config["output-dir"] + "/step{}_solvation.gro".format(current_step - 1), config["output-dir"] + "/_ion.gro")

        # Generate ions
        for ion_pair in ions:

            positive = ion_pair["positive-ion"]
            negative = ion_pair["negative-ion"]
            concentration = ion_pair["concentration-M"]
            neutralise = ion_pair["neutralise"]

            execute_gmx("grompp -f $OUT_DIR$/arbitrary.mdp -c $OUT_DIR$/_ion.gro -o $OUT_DIR$/_ionisation.tpr -p $OUT_DIR$/topol.top", step=current_step)

            precommand = "printf '%s\n' 13"

            os.remove(config["output-dir"] + "/_ion.gro")

            if neutralise:
                execute_gmx("genion -s $OUT_DIR$/_ionisation.tpr -p $OUT_DIR$/topol.top -pname {} -nname {} -o $OUT_DIR$/_ion.gro -conc {} -noneutral".format(positive, negative, str(concentration)), step=current_step, precommand=precommand)

            else:
                execute_gmx("genion -s $OUT_DIR$/_ionisation.tpr -p $OUT_DIR$/topol.top -pname {} -nname {} -o $OUT_DIR$/_ion.gro -conc {} -neutral".format(positive, negative, str(concentration)), step=current_step, precommand=precommand)

            os.remove(config["output-dir"] + "/_ionisation.tpr")

        os.rename(config["output-dir"] + "/_ion.gro", config["output-dir"] + "/step{}_ions.gro".format(current_step))

        last_out = "step{}_ions.gro".format(current_step)

logger.log(2, "Generating a universal index file")
execute_gmx("make_ndx -f $OUT_DIR$/{} -o $OUT_DIR$/index.ndx".format(last_out), step=current_step, precommand="printf '%s\n' q")

# Perform steps
logger.log(2, "Preparing runs - total of {} steps".format(len(config["steps"])))

for run_step in config["steps"]:

    reference = run_step["reference"]
    mdp_shortname = reference.replace(".mdp", "")

    # Copy over the MDP file and edit it
    mdp_data = mdpreader.mdp_to_json(open("templates/" + reference).read())

    # Edit configuration
    if "dt" not in mdp_data.keys():
        mdp_data["dt"] = None

    # In picoseconds
    mdp_data["dt"] = run_step["femtoseconds-per-step"] / 1000

    if "nsteps" not in mdp_data.keys():
        mdp_data["nsteps"] = None

    simulation_length = run_step["simulation-length-picoseconds"]
    total_iterations = math.ceil(simulation_length / mdp_data["dt"])

    mdp_data["nsteps"] = total_iterations

    # Other edits
    if run_step["additional-configs"] != None:

        # continue
        for addable in run_step["additional-configs"].keys():

            mdp_data[addable] = run_step["additional-configs"][addable]

    out_mdp = mdpreader.json_to_mdp(mdp_data)

    mdp_file = "step{}_".format(current_step) + reference
    open(output_directory + "/" + mdp_file, "w+").write(out_mdp)

    execute_gmx("grompp -f $OUT_DIR$/{} -c $OUT_DIR$/{} -n $OUT_DIR$/index.ndx -p $OUT_DIR$/topol.top -o $OUT_DIR$/step$STEP$_{} -maxwarn 10".format(mdp_file, last_out, mdp_shortname), step=current_step)

    simulation_name = "step" + str(current_step) + "_" + mdp_shortname

    logger.log(2, "Performing simulation {}, total {} iterations".format(simulation_name, total_iterations))

    # Perform the MD
    execute_gmx("mdrun -deffnm $OUT_DIR$/{}".format(simulation_name), step=current_step)

    logger.log(2, "Simulation {} finished".format(simulation_name))

    # Simulation output

    last_out = simulation_name + ".gro"
    current_step += 1
