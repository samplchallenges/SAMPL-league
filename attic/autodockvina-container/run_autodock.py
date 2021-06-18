import os
import click
import docker


# taken from https://hub.docker.com/r/taccsciapps/autodock-vina
# found using "docker run taccsciapps/autodock-vina -e"


def create_sh_file(config_file, pdbqt_in, pdbqt_out):
	''' creates a file called "dock.sh" which has the run command to run autodock fina 
		using taccsciapp/autodock-vina container
	'''
	ofile = open("./dock.sh", 'w')
	ofile.write("#!/bin/bash\n")
	ofile.write("mkdir -p ./out\n")
	runcmd = "/usr/bin/time -p vina --config {} --ligand {} --out ./out/{} > /dev/null;".format(config_file, pdbqt_in, pdbqt_out)
	ofile.write(runcmd + "\n")
	ofile.write("score=$(/bin/grep -m 1 \"REMARK VINA RESULT:\" ./out/ZINC00000567.pdbqt | /usr/bin/awk '{ print $4 }');\n")
	ofile.write("echo \"./out/ZINC00000567.pdbqt $score\" >> ./out/tempResults\n")
	ofile.close()
	os.system("cat ./dock.sh")
	print("done")


@click.command()
@click.argument("config_file")
@click.argument("pdbqt_in")
@click.argument("pdbqt_out")
def autodock(config_file, pdbqt_in, pdbqt_out):
	#create_sh_file(config_file, pdbqt_in, pdbqt_out)

	#os.system("docker run -it --rm -v $(pwd):/data taccsciapps/autodock-vina /dock.sh")
	#client = docker.from_env()
	#client.containers.run("taccsciapps/autodock-vina", "/examples/demo.sh")


	#os.system("/usr/bin/time -p vina")
	#os.system("./vina --config {} --ligand {} --out {} --score-only".format(
	#	config_file, pdbqt_in, pdbqt_out))

	#os.system("./vina")

if __name__ == "__main__":
	print("starting")
	os.system("docker run -it --rm -v $(pwd):/data taccsciapps/autodock-vina /dock.sh")

