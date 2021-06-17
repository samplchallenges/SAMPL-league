import os
import click
import docker

from docker.types import Mount

DEBUG = True


# taken from https://hub.docker.com/r/taccsciapps/autodock-vina
# found using "docker run taccsciapps/autodock-vina -e"



def make_sh_file(cwd, shfile, config_file, pdbqt_in, pdbqt_out):
	''' creates a file called "dock.sh" which has the run command to run autodock fina 
		using taccsciapp/autodock-vina container
	'''
	outdir = cwd + "/out/"
	
	config_path = "/examples/" + config_file
	pdbqtin_path = "/examples/" + pdbqt_in
	pdbqtout_path = outdir + pdbqt_out

	ofile = open(shfile, 'w')
	ofile.write("#!/bin/bash\n")
	ofile.write("mkdir -p {}\n".format(outdir))
	runcmd = "/usr/bin/time -p vina --config {} --ligand {} --out {};".format(config_path, pdbqtin_path, pdbqtout_path)
	ofile.write(runcmd + "\n")
	#ofile.write("score=$(/bin/grep -m 1 \"REMARK VINA RESULT:\" ./out/ZINC00000567.pdbqt | /usr/bin/awk '{ print $4 }');\n")
	#ofile.write("echo \"./out/ZINC00000567.pdbqt $score\" >> ./out/tempResults\n")
	ofile.close()
	if DEBUG: 
		os.system("cat {}/dock.sh".format(cwd))
		print()

	


@click.command()
@click.argument("config_file") # strings whether or not input as "config.in" or config.in on command line
@click.argument("pdbqt_in")
@click.argument("pdbqt_out")
def autodock(config_file, pdbqt_in, pdbqt_out):	

	if DEBUG:
		print("CONFIG_INPUT:   type: {}   input: {}".format(type(config_file), config_file))
		print("PDBQTI_INPUT:   type: {}   input: {}".format(type(pdbqt_in), pdbqt_in))
		print("PDBQTO_INPUT:   type: {}   input: {}".format(type(pdbqt_out), pdbqt_out))
		print()
	cwd = os.getcwd()
	shfile = "examples/dock.sh"
	make_sh_file(cwd, shfile, config_file, pdbqt_in, pdbqt_out)
	os.system("chmod +x {}".format(shfile))
	print("chmodded\n")
	client = docker.from_env()
	output = client.containers.run("taccsciapps/autodock-vina", "dock.sh")
	print(output.decode("utf-8"))


	#os.system("/usr/bin/time -p vina")
	#os.system("./vina --config {} --ligand {} --out {} --score-only".format(
	#	config_file, pdbqt_in, pdbqt_out))

	#os.system("./vina")

if __name__ == "__main__":
	#autodock()
	print("program started")
	cwd = os.getcwd()
	mt = Mount("taccsciapps/autodock-vina", cwd)
	client = docker.from_env()

	container = client.api.create_container('taccsciapps/autodock-vina', host_config=client.api.create_host_config(binds=["/Users/megosato/SAMPL-league/attic/autodockvina-container/ad-test:/data"]))
	
	print("running container id:", container)
	# client.containers.run(container.get('Id'))     			does not work
	# client.containers.run(container, "./dock.sh")  			does not work
	# client.containers.run(container_id["Id"], "./dock.sh")	does not work
	# print(output)
