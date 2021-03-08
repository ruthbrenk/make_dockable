#!/usr/bin/env python

#load smiles in taut_smiles table 

import string,os, sys, my_mysql3 as mysql, datetime, time
import pymysql

if len(sys.argv) != 6 :
	print ('mysql_load_taut_states.py  <username> <password> <smiles.smi> <project> <initial> >')
	sys.exit(1)

username = sys.argv[1]
password = sys.argv[2]
data_file = open(sys.argv[3], 'r')
project = sys.argv[4]
initial = sys.argv[5]

#date_of_update 
utc_datetime = datetime.datetime.utcnow()
date_of_update  = str(utc_datetime.strftime("%Y-%m-%d"))
#-------------------------------------
#-------------------------------------

  
conn=mysql.connect2server(password,username, 'purchasable')  			  
			  
cursor = conn.cursor ()

#get table names

command = 'select tautomerized_smiles,protonated_smiles from projects where project = "' + project + '"'
print (command)
cursor.execute(command)
results=cursor.fetchall()
if len(results) == 0:
	print ('project does not exist')
	sys.exit(1)

taut_table = results[0][0]
prot_table = results[0][1]
print (taut_table, prot_table)

#get comp_id, smiles_field from table
command = "show fields from " + taut_table
cursor.execute(command)
results=cursor.fetchall()
taut_id_field = results[0][0]
id_field = results[1][0]
smiles_field = results[2][0]
prot_id_field = results[3][0]

print (taut_id_field, smiles_field, id_field, prot_id_field)


print ("insert data")

counter = 0

new_smiles = {}

old_id = ''

for i in data_file.readlines():
	#smi, identifier, state = string.split(i.strip(),'\t')
	smi, prot_id = i.strip().split('\t')
	smi = mysql.encode_smiles(smi) 

	#get id
	command = "select " + id_field + " from " + prot_table + " where " + prot_id_field + "= " + prot_id
	#print command
	cursor.execute(command)
	results = cursor.fetchall()
	id = str(results[0][0])

	#delete old entries
	if old_id != prot_id:
	    	command = "delete from " + taut_table + " where " + prot_id_field + " =" + prot_id
	    	#print command
	    	cursor.execute(command)

	#insert into tautomerized_smiles
	command = "INSERT INTO " + taut_table + " (" + id_field + ", " + prot_id_field + ", " + smiles_field + ", initial, date_of_update) VALUES (" + id + ", " + prot_id + ", '" + smi + "', '" + initial + "', '"  + date_of_update + "')"
	#pint command

	error =True

	while error: #catch Deadlock

		try:
			cursor.execute(command)
			error = False

		# NB : you won't get an IntegrityError when reading
		except (pymysql.Error, pymysql.Warning) as e:
			if e[0] == 1213:
				print (e[0], e[1])
				time.sleep(15)
		
 
	old_id = prot_id

cursor.close ()
#next line necessary for InnoDB tables, otherwise updates don't show up
conn.commit()
conn.close ()
