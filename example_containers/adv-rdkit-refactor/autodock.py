from rdkit import Chem
from rdkit.Chem import AllChem
import os

class Autodock():
	def __init__(self, python_path: str, utility_path: str, vina_path: str, receptor: str, box: tuple):
		self.python_path = python_path 		
		self.utility_path = utility_path	# should end without '/'
		self.vina_path = vina_path

		self.receptor_file_path = receptor
		self.sz_x, self.sz_y, self.sz_z, self.c_x, self.c_y, self.c_z = box


	def __str__(self):
		return f"Autodock: Receptor=\"{self.receptor_file_path}\", Box=({self.sz_x},{self.sz_y},{self.sz_z},{self.c_x},{self.c_y},{self.c_z}), Python=\"{self.python_path}\", Utilities=\"{self.utility_path}\", Vina=\"{self.vina_path}\""


	def _get_utility_cmd(self, pyfile: str):
		return f"{self.python_path} {self.utility_path}/{pyfile}"

	def _make_config_file(self, exhaustiveness, num_modes, flex, config_path):
		''' makes the configuration file which is an input to the vina 
		    command
		'''
		configfile = open(config_path, "w")
		configfile.write(f"receptor = {self.recprep_pdbqt}\n")

		configfile.write(f"center_x = {self.c_x}\n")
		configfile.write(f"center_y = {self.c_y}\n")
		configfile.write(f"center_z = {self.c_z}\n")

		configfile.write(f"size_x = {self.sz_x}\n")
		configfile.write(f"size_y = {self.sz_y}\n")
		configfile.write(f"size_z = {self.sz_z}\n")

		if exhaustiveness != None:
			configfile.write(f"exhaustiveness = {exhaustiveness}\n")
		if num_modes != None:
			configfile.write(f"num_modes = {num_modes}")
		
		configfile.close()


	def _save_highest_score(self, ligdock_path: str, lighighscore_path: str):
		''' writes the highest scoring docked pose into its own file
		'''

		infile = open(ligdock_path, 'r')
		outfile = open(lighighscore_path, 'w')
		
		for line in infile.readlines():
			if "MODEL" in line:
				continue
			if "ENDMDL" in line:
				break
			outfile.write(line)
				
		outfile.close()
		infile.close()

	def _get_score(self, score_path):
		with open(score_path) as scoref:
			for line in scoref:
				if "Aff" in line:
					return float(line.split()[1])

	@staticmethod
	def charge_ligand(smiles, outfile: str):
		mol = Chem.MolFromSmiles(smiles)
		mol2 = Chem.AddHs(mol)
		AllChem.EmbedMolecule(mol2)

		# save the charged molecule as an sdf
		ostream = Chem.SDWriter(outfile)
		ostream.write(mol2)

	@staticmethod
	def sdf_to_pdbqt(sdf_file: str, pdbqt_file: str):
		os.system(f"obabel {sdf_file} -O {pdbqt_file} 2>/dev/null 1>/dev/null")

	def pdbqt_to_pdb(self, pdbqt_path, pdb_path):
		cmd = self._get_utility_cmd("pdbqt_to_pdb.py")
		#os.system(f"cut -c-66 {pdbqt_path} > {pdb_path}")
		os.system(f"{cmd} -f {pdbqt_path} -o {pdb_path}")

	def prep_ligand(self, ligchg_pdbqt: str, ligprep_pdbqt: str):
		cmd = self._get_utility_cmd("prepare_ligand4.py")
		os.system(f"{cmd} -l {ligchg_pdbqt} -o {ligprep_pdbqt}")


	def prep_receptor(self, recprep_pdbqt: str):
		cmd = self._get_utility_cmd("prepare_receptor4.py")
		self.recprep_pdbqt = recprep_pdbqt
		os.system(f"{cmd} -r {self.receptor_file_path} -o {recprep_pdbqt} > /dev/null")


	def dock(self, ligprep_pdbqt: str, ligdock_pdbqt: str, lighighscore_pdbqt: str, lighighscore_pdb: str, exhaustiveness: int, num_modes: int, flex):
		config_path = "/tmp/config.txt"
		score_path = "/tmp/score.txt"

		self._make_config_file(exhaustiveness, num_modes, flex, config_path)
		os.system(f"{self.vina_path} --config {config_path} --ligand {ligprep_pdbqt} --out {ligdock_pdbqt} > /dev/null")

		self._save_highest_score(ligdock_pdbqt, lighighscore_pdbqt)
		self.pdbqt_to_pdb(lighighscore_pdbqt, lighighscore_pdb)

		os.system(f"{self.vina_path} --config {config_path} --ligand {lighighscore_pdbqt} --out /tmp/rescore.pdbqt --score_only | grep 'Affinity:' > {score_path}")
		score = self._get_score(score_path)
		return score

	def save_receptor_pdb(recprep_pdb: str):
		Autodock.pdbqt_to_pdb(recprep_pdbqt, recprep_pdb)



