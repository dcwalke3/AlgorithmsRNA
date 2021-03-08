'''
By: Dakota Walker
Class: Algorithms
Assignment: Protein Translation
'''

#Gets file containing the Raw DNA code.
with open('RawProtein.txt', 'r') as file:
    #Makes the multiple line text file into one long string to be split later.
    data = file.read().replace('\n', '')
    

protienCountGlobal = {"U":0, "A":0, "G":0, "C":0, "Total":0}

#Translates the data to RNA.
Translate = {"T":"A", "A":"U", "G":"C", "C":"G"}
newData = ""
for c in range(len(data)):
    newData += Translate[data[c]]

# Seperates things into a group of 3.
n=3
seperatedData = [newData[i:i+n] for i in range(0, len(newData), n)]


#Used to count protein Frequency per protein.
proteinFrequencyCounter = {"A":0, "C":0, "G":0, "U":0, "Total":0}

#Counts the number of certain types of protein
proteinStrandFrequencyGlobal = {"Asp":0, "Asn":0, "Arg":0, "Tyr":0, "Gln":0,
                                "Trp":0, "Glu":0, "His":0, "Cys":0, "Pro":0,
                                "Ile":0, "STOP":0, "Thr":0, "Met":0, "Leu":0,
                                "Ala":0, "Gly":0, "Val":0, "Ser":0, "Phe":0,
                                "Lys":0 ,"Total":0}

#Each number contains a protein Frequency Dictionary, but it is done at
NumbersCounted = {}


#Conversion key for the proteins to help me
ProteinConversionKey = {"UUU":"Phe", "UUC":"Phe", "UUG":"Leu", "UUA":"Leu",
                        "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser",
                        "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
                        "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
                        "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",
                        "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
                        "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
                        "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
                        "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
                        "ACU":"Thr", "ACC":"Thr", "ACG":"Thr", "ACA":"Thr",
                        "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
                        "AGU":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
                        "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
                        "GCU":"Ala", "GCA":"Ala", "GCC":"Ala", "GCG":"Ala",
                        "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
                        "GGU":"Gly", "GGC":"Gly", "GGG":"Gly", "GGA":"Gly",}


translatedData = []
tempList=[]

#Made this just to make counting ACGUs easier.
tempRawString=""
dictIndex=1
for i in range(len(seperatedData)): 
    if (ProteinConversionKey[seperatedData[i]]=="STOP"):
        
        proteinStrandFrequencyGlobal[ProteinConversionKey[seperatedData[i]]]+=1
        proteinStrandFrequencyGlobal["Total"]+=1
        
        tempList.append(ProteinConversionKey[seperatedData[i]])
        translatedData.append(tempList)
        tempList=[]

        #Used to count the number of ACGUs with each protein strand.
        NumbersCounted[dictIndex] = {"A":0, "U":0, "G":0, "C":0, "Total":0}
        tempRawString+=seperatedData[i]
        for z in range(len(tempRawString)):
            NumbersCounted[dictIndex][tempRawString[z]]+=1
            NumbersCounted[dictIndex]["Total"]+=1    
        dictIndex+=1
        tempRawString=""
    else:
        proteinStrandFrequencyGlobal[ProteinConversionKey[seperatedData[i]]]+=1
        proteinStrandFrequencyGlobal["Total"]+=1
        tempList.append(ProteinConversionKey[seperatedData[i]]+"-")
        tempRawString+=seperatedData[i]

#Returns all protein strands into strings instead of individual list elements.            
stringData=""
for i in range(len(translatedData)):
    stringData = "".join(translatedData[i])
    translatedData[i]=stringData

'''
Writes information to a new .txt file (if file does not exist else will rewrite over).
Had to use a lot of format.
If wanted you can remove ProteinTranslation.txt and it will still work.
'''
results = open("ProteinTranslation.txt", "w")

results.writelines([
    "Initial Data: "+data+"\n",
    "Translated Data: "+newData+"\n\n",
])

for i in range(len(NumbersCounted)):
    results.writelines([
        translatedData[i],
    "\nA  {:.2%}\n".format((NumbersCounted[i+1]["A"]/NumbersCounted[i+1]["Total"])),
    "C  {:.2%}\n".format((NumbersCounted[i+1]["C"]/NumbersCounted[i+1]["Total"])),
    "G  {:.2%}\n".format((NumbersCounted[i+1]["G"]/NumbersCounted[i+1]["Total"])),
    "U  {:.2%}".format((NumbersCounted[i+1]["U"]/NumbersCounted[i+1]["Total"]))+"\n\n",
    ])
results.writelines([
    "----Protein Frequencies----\n",
])

for key in proteinStrandFrequencyGlobal:
    if(proteinStrandFrequencyGlobal[key]>0):    
        results.writelines([
        ""+key+"   {:.2%}\n".format((proteinStrandFrequencyGlobal[key]/proteinStrandFrequencyGlobal["Total"])),
        ])
results.writelines([ 
    "\n\n",
    "Total Protiens: "+str(len(translatedData))
])
results.close()