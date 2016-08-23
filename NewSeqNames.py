from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests, os, re
from openpyxl import load_workbook

fastafile = os.path.join(os.path.abspath("../.."), u"Python_fasta_files/dloop/")
wb = load_workbook(filename = "/Users/katalinabobowik/Documents/MEEP_Lab/Data_Curation/Data_Curation_dloop_V2.xlsx")


# regex to find non-numeric characters
regexp = re.compile("[^0-9]")


for name in wb.get_sheet_names():

    count = 0
    sheet = list(wb.get_sheet_by_name(name))

    header = sheet.pop(0)
    header = map((lambda x: x.value), header)

    try:
        accession = header.index("Accession Number")
        radiocarbon = header.index("Mean Age")
        sigma = header.index("Sigma")

    except:
        print "{} is missing accession number or mean age".format(name)
        print header

    new_file_name = fastafile + name
    records = list(SeqIO.parse(new_file_name + ".fasta", "fasta"))
    fout = open("{}.fa".format(name), "w")

    new_records = []

    for row in sheet:

        try:
            if not regexp.search(str(row[radiocarbon].value)) and row[accession].value in records[count].id:
                #print "Writing"
                #print records[count].id

                try:
                    #records[count].id = str(row[accession].value) + "_" + str(row[radiocarbon].value)
                    #write_me = str(row[accession].value) + "_" + str(row[radiocarbon].value) + " " + str(records[count].seq)

                    new_id = str(row[accession].value) + "_" + str(row[radiocarbon].value) + "_" + str(row[sigma].value)
                    new_record = SeqRecord(records[count].seq, new_id, '', '')

                    new_records.append(new_record)

                except:
                    print "This is an except error"
                    print name
                    print records[count].id
                    print row[accession].value, row[radiocarbon].value

            else:

                """
                print "Else error"
                print name
                print row[accession].value
                print records[count].id
                """
        except:

            print "This is the outside except error"
            print name
            print records[count].id
            print row[accession].value, row[radiocarbon].value



        count+=1

    SeqIO.write(new_records, fout, 'fasta')

    fout.close()

