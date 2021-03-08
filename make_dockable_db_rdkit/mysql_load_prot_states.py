#!/usr/bin/env python

#load smiles in prot_smiles table 

import string,os, sys, my_mysql3 as mysql, datetime, time
import pymysql

if len(sys.argv) != 6 :
	print ('mysql_load_prot_states.py  <username> <password> <smiles.smi> <project> <initial> >')
	sys.exit(1)

username = sys.argv[1]
password = sys.argv[2]
data_file = open(sys.argv[3], 'r')
project = sys.argv[4]
initial = sys.argv[5]

#date_of_update 
utc_datetime = datetime.datetime.utcnow()
date_of_update  = str(utc_datetime.strftime("%Y-%m-%d"))

counter = 0

#-------------------------------------
#-------------------------------------

  
conn=mysql.connect2server(password,username,'purchasable')  			  
			  
cursor = conn.cursor ()



#get table names

command = 'select protonated_smiles from projects where project = "' + project + '"'
print (command)
cursor.execute(command)
results=cursor.fetchall()
if len(results) == 0:
	print ('project does not exist')
	sys.exit(1)

prot_table = results[0][0]
print (prot_table)

#get comp_id, smiles_field from table
command = "show fields from " + prot_table
cursor.execute(command)
results=cursor.fetchall()
prot_id_field = results[0][0]
smiles_field = results[2][0]
id_field = results[1][0]

print (prot_id_field, smiles_field, id_field)


print ("insert data")


new_smiles = {}

old_id = ''

for i in data_file.readlines():

    smi, identifier = i.strip().split('\t')

    if identifier != old_id:
    	command = "delete from " + prot_table + " where " + id_field + "=" + identifier 
    	#print (command)
    	cursor.execute(command)   
    old_id = identifier

    smi = mysql.encode_smiles(smi) 	

    command = "INSERT INTO " + prot_table + " (" + id_field + ", " + smiles_field + ", initial, date_of_update) VALUES (" + identifier + ", '" + smi  + "', '" + initial + "', '"  + date_of_update + "')"
    #print command

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

cursor.close ()
#next line necessary for InnoDB tables, otherwise updates don't show up
conn.commit()
conn.close ()
