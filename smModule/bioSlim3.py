import MySQLdb as mdb
import csv
import sys
import smModule.misc as misc

with open("configNew.txt") as c:
	config = misc.readConfig(c.readlines())



def flatl(l):
	"""
	Flattening lists made easy!
	Or unpacking tuples made easy!
	Try:
	a = ['a','b',[1,2,3,4]]
	print flatl(a)
	"""
	return([item for subl in l for item in subl])

def merge_intv(intervals, distance = 6):
	"""
	An algorithm to sort tuples with (start, end) coordinates
	the last element is compared with the upper element and a new one is generated
	on the fly in case the coordinates are adjacent to each other
	Default argument for distance is 6, can be modified by providing a dist argument
	"""
	sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
	merged = []
	global lower

	for higher in sorted_by_lower_bound:
		if not merged:
			merged.append(higher)
		else:
			lower = merged[-1]
			# test for intersection between lower and higher:
			# we know via sorting that lower[0] <= higher[0]
		if all([lower[1]-lower[0] < 300, higher[1]-higher[0] < 300]):
				if abs(higher[0]-lower[1]) < distance: #higher[0] <= lower[1]:
					upper_bound = max(lower[1], higher[1])
					merged[-1] = (lower[0], upper_bound)  # replace by merged interval
				elif abs(higher[0]-lower[0]) < distance: #higher[0] <= lower[1]:
					upper_bound = max(lower[1], higher[1])
					merged[-1] = (lower[0], upper_bound)
				elif higher[1]-lower[1] < 0: #merging elements already included (due to different algorithms)
					upper_bound = max(lower[1], higher[1])
					merged[-1] = (lower[0], upper_bound)
				else:
					merged.append(higher)
		else:
			merged.append(higher)
	return merged


def read_fasta(ff):
	"""
	Actually from biopython
	seqDict = {}
	with open('/home/seth/domainPhylogeny/fastaFiles/hybridsOnlypubData.fasta') as ff:
	for name, seq in read_fasta(ff):
	seqDict[name] = seq
	"""
	name, seq = None, []
	for line in ff:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line[1:], []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


def cleanProtSeq(seq):
	if '*' in seq:
		seq = seq.replace('*','x')
	if seq.endswith('x'):
		seq = seq[:-1]

	return seq

def tupleToFasta(db, fileName, maxChars = 80):
	"""
	Provide a simple dictionary with identifier as key and sequence as value to write it to
	a fasta file
	"""
	outfile = open(fileName, 'w')
	counter = 1

	for (key, seq) in db:
		seq = cleanProtSeq(seq)
		outfile.write('>' + str(key) +'\n')
		counter+=1
		for i in range(0, len(seq), maxChars):
			outfile.write(seq[i:i+maxChars]+"\n")
	outfile.close()

def dictToFasta(db, fileName, maxChars = 80):
	"""
	Provide a simple dictionary with identifier as key and sequence as value to write it to
	a fasta file
	"""
	outfile = open(fileName, 'w')
	counter = 1

	for key, seq in db.items():
		outfile.write('>' + str(key) +'\n')
		counter+=1
		for i in range(0, len(seq), maxChars):
			outfile.write(seq[i:i+maxChars]+"\n")
	outfile.close()





def dbFetch(query):
	db = mdb.connect(host=config['host'], user=config['user'], passwd=config['passwd'], db=config['db'])
	cursor = db.cursor()
	cursor.execute(query)
	data = list(cursor.fetchall())
	return(data)
# Why doesnt this work?

def dbwHeader(query):
	db = mdb.connect(host=config['host'], user=config['user'], passwd=config['passwd'], db=config['db'])
	cursor = db.cursor()
	cursor.execute(query)
	num_fields = len(cursor.description)
	field_names = tuple([i[0] for i in cursor.description])
	print("Mysql query yielded info on")
	print(field_names)
	data = list(cursor.fetchall())
	print("Number of rows:")
	print(len(data))
	return([field_names, data])


def iprFileReader(iprFile):
	print("Reading interpro file")
	npDomains = {}
	with open(iprFile) as tsv:
		handle = csv.reader(tsv, delimiter = '\t')
		for line in handle:
			if len(line) < 13: # exlcuding annotations without interpro entry
				continue
			name, md5, unknown, domainDb, pf_id, pf_desc, start, end, score, unknown2, date, ipr_id, ipr_desc = line
			if name not in npDomains.keys():
				npDomains[name] = []
			npDomains[name].append((int(start),int(end)))
	return(npDomains)



def proteinDl(combinedId):
    """
    Downloads all secondary metabolite proteins (the ones like, e.g. PKS, NRPS, etc.) for a given organism.
    """
    print("Downloading secondary metabolite proteins")

    proteins = bio.dbFetch("""
    SELECT torg.name, torg.org_id, proteins.prot_seqkey, sp.sm_short, proteins.prot_seq FROM (SELECT * FROM organism WHERE name IN ('%s')) torg
    JOIN smurf_papa AS sp ON torg.org_id = sp.org_id AND sp.sm_short != 'none'
    JOIN proteins ON sp.org_id = proteins.org_id AND sp.sm_protein_id = proteins.prot_seqkey;
    """ % "','".join(orgs) )

    proteins = [(org, org_id, protein_id, sm_short, bio.cleanProtSeq(seq.decode("UTF-8"))) for org, org_id, protein_id, sm_short, seq in proteins]

    return(proteins)


class smProt:
	'Class for secondary metabolite associated proteins'

	def __init__(self, name):
	  self.name = name

	def setSeq(self, seq):
		try:
			isinstance(seq, str)
			self.seq = seq
		except Exception:
			print("Error, sequence must be string")

	def setDomains(self, domains, minSize = 100, distance = 100):
		try:
			print(isinstance(domains, list))
			isinstance(domains, list)
			domains = bio.merge_intv(domains, distance = distance) # Joining domains together
			domains = [i for i in domains if i[1]-i[0] > minSize] # Only taking domains over 250
			self.domains = domains
		except Exception:
			print ("Error: Domains must be formatted as list")

	def fullSeq(self):
		return self.seq

	def getDomains(self):
		return self.domains

	def getDomainList(self):
		return [[self.name] + list(d) for d in self.domains]

	def getDomSeqs(self):
		l = [ ]
		print(self.seq)
		for start, end in self.domains:
			dseq = self.seq[start-1:end]
			l.append(dseq)
		return(l)

	def displayProtein(self):
		return "Name : ", self.name,  ", Domains: ", self.domains, ", Sequence Length: ", len(self.seq)
