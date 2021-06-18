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
	outdir = cwd 


	ofile = open(shfile, 'w')
	ofile.write("#!/bin/bash\n")
	
	runcmd = "/usr/bin/time -p vina --config {} --ligand {} --out {};".format(config_file, pdbqt_in, pdbqt_out)
	ofile.write(runcmd + "\n")

	ofile.close()


	if DEBUG: 
		os.system("cat {}".format(shfile))
		print()

	
def run_autodock(shfile):
	cwd = os.getcwd()
	client = docker.from_env()
	client.containers.run('taccsciapps/autodock-vina', "./dock.sh", volumes={cwd: {'bind': '/data', 'mode': 'rw'}})



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
	shfile = "dock.sh"
	make_sh_file(cwd, shfile, config_file, pdbqt_in, pdbqt_out)

	outdir = cwd
	#print("making directory", outdir)
	#os.mkdir(outdir)


	os.system("chmod +x {}".format(shfile))
	os.system("docker run -v $(pwd):/data taccsciapps/autodock-vina {}".format(shfile))

	
	run_autodock(shfile)

if __name__ == "__main__":
	print("program started")
	# CLI Run Command: "docker run -it --rm -v $(pwd):/data taccsciapps/autodock-vina dock.sh"
	# os.system("docker run -it --rm -v $(pwd):/data taccsciapps/autodock-vina dock.sh")



	autodock()
	



	
	#container = client.api.create_container('taccsciapps/autodock-vina', host_config=client.api.create_host_config(binds=["{}:/data".format(cwd)]))
	#print(type(container))
	






	#cwd = os.getcwd()
	#mt = Mount("taccsciapps/autodock-vina", cwd)
	#client = docker.from_env()

	#container = client.api.create_container('taccsciapps/autodock-vina', host_config=client.api.create_host_config(binds=["/Users/megosato/SAMPL-league/attic/autodockvina-container/ad-test:/data"]))
	
	#print("running container id:", container)
	# client.containers.run(container.get('Id'))     			does not work
	# client.containers.run(container, "./dock.sh")  			does not work
	# client.containers.run(container_id["Id"], "./dock.sh")	does not work
	# print(output)
