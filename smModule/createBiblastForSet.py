import MySQLdb as mdb
import misc
import sys


with open("../configNew.txt") as c:
	config = misc.readConfig(c.readlines())

with open("../nigri_set.txt") as orgf:
	orgSet = [org.strip() for org in orgf.readlines()]
# orgsForQuery(orgSet)

def orgsForQuery(orgs):
	return("('"+("','").join(orgs)+"')")

# usually would be blastdb
def biblastMysqlHandler(cursor, smtable,orgSet):
	query = """CREATE TABLE %s AS
	SELECT * FROM blastdb.biblast WHERE q_org IN %s AND h_org IN %s;""" % (smtable, orgsForQuery(orgSet), orgsForQuery(orgSet))

	print("Executing query")
	print(query)
	try:
		cursor.execute(query)
	except Exception as e:
		raise print("Cannot create table")
		db.close()
	print("Done")


	indexQueries = ["CREATE INDEX i_qorg ON %s (q_org);" % smtable,
	"CREATE INDEX i_horg ON %s (h_org);" % smtable,
	"CREATE INDEX i_q_seqkey ON %s (q_seqkey);"% smtable,
	"CREATE INDEX i_h_seqkey ON %s (h_seqkey);"% smtable]

	print("Creating indexes")
	for query in indexQueries:
		try:
			cursor.execute(query)
		except Exception as e:
			raise print("Cannot run indexes")
			db.close()

	return("Finished running query and indexing")



def createBidirSet(smtable = "null", orgSet = "null"):

	if any([smtable == "null", orgSet == "null"]):
		raise ValueError("Please make sure to provide a table name and an organism set")

	db = mdb.connect(host=config['host'], user=config['user'], passwd=config['passwd'], db=config['db'])
	cursor = db.cursor()

	query = "SHOW TABLES"
	cursor.execute(query)
	tables = [item[0] for item in  cursor.fetchall()]

	if smtable in tables:
		print("Table %s already exists. If you performed a run on the same set before, keep it.\
		If it is on another set or you don't know, delete it.\n" % smtable)

		smurfTempAns = input("Delete table %s and generate new one? (y/n)\n" % smtable).lower()

		if smurfTempAns.startswith('y'):
			print("Deleting table %s" % smtable)
			query = """DROP TABLE %s;""" % smtable
			cursor.execute(query)

			biblastMysqlHandler(cursor,smtable, orgSet)

			print("Done")
		else:
			print("Keeping %s \n" % smtable)
	else:
		print("Creating new table %s" % smtable)

		biblastMysqlHandler(cursor,smtable, orgSet)

	print("Done")
	db.close()

if __name__ == '__main__':
	createBidirSet(smtable = "publication_nigri_biblast", orgSet = orgSet)
