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
newData = ""
for c in range(len(data)):
    char = ""
    if data[c] == "T":
        char = "A"
        protienCountGlobal["A"]+=1
    elif data[c] == "A":
        char = "U"
        protienCountGlobal["U"]+=1
    elif data[c] == "G":
        char = "C"
        protienCountGlobal["C"]+=1
    elif data[c] == "C":
        char = "G"
        protienCountGlobal["G"]+=1
    newData += char
    protienCountGlobal["Total"]+=1


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
ProteinConversionKey = {"UUU":"Phe", "UUC":"Phe", "UUU":"Leu", "UUU":"Leu",
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
    
    if(seperatedData[i] in ProteinConversionKey):
        
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
                if(tempRawString[z]=="A"):
                    NumbersCounted[dictIndex]["A"]+=1
                elif(tempRawString[z]=="U"):
                    NumbersCounted[dictIndex]["U"]+=1
                if(tempRawString[z]=="G"):
                    NumbersCounted[dictIndex]["G"]+=1
                if(tempRawString[z]=="C"):
                    NumbersCounted[dictIndex]["C"]+=1
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
    "Asp   {:.2%}\n".format((proteinStrandFrequencyGlobal["Asp"]/proteinStrandFrequencyGlobal["Total"])),
    "Asn   {:.2%}\n".format((proteinStrandFrequencyGlobal["Asn"]/proteinStrandFrequencyGlobal["Total"])),
    "Arg   {:.2%}\n".format((proteinStrandFrequencyGlobal["Arg"]/proteinStrandFrequencyGlobal["Total"])),
    "Tyr   {:.2%}\n".format((proteinStrandFrequencyGlobal["Tyr"]/proteinStrandFrequencyGlobal["Total"])),
    "Gln   {:.2%}\n".format((proteinStrandFrequencyGlobal["Gln"]/proteinStrandFrequencyGlobal["Total"])),
    "Glu   {:.2%}\n".format((proteinStrandFrequencyGlobal["Glu"]/proteinStrandFrequencyGlobal["Total"])),
    "Trp   {:.2%}\n".format((proteinStrandFrequencyGlobal["Trp"]/proteinStrandFrequencyGlobal["Total"])),
    "His   {:.2%}\n".format((proteinStrandFrequencyGlobal["His"]/proteinStrandFrequencyGlobal["Total"])),
    "Cys   {:.2%}\n".format((proteinStrandFrequencyGlobal["Cys"]/proteinStrandFrequencyGlobal["Total"])),
    "Pro   {:.2%}\n".format((proteinStrandFrequencyGlobal["Pro"]/proteinStrandFrequencyGlobal["Total"])),
    "Ile   {:.2%}\n".format((proteinStrandFrequencyGlobal["Ile"]/proteinStrandFrequencyGlobal["Total"])),
    "STOP  {:.2%}\n".format((proteinStrandFrequencyGlobal["STOP"]/proteinStrandFrequencyGlobal["Total"])),
    "Thr   {:.2%}\n".format((proteinStrandFrequencyGlobal["Thr"]/proteinStrandFrequencyGlobal["Total"])),
    "Met   {:.2%}\n".format((proteinStrandFrequencyGlobal["Met"]/proteinStrandFrequencyGlobal["Total"])),
    "Leu   {:.2%}\n".format((proteinStrandFrequencyGlobal["Leu"]/proteinStrandFrequencyGlobal["Total"])),
    "Ala   {:.2%}\n".format((proteinStrandFrequencyGlobal["Ala"]/proteinStrandFrequencyGlobal["Total"])),
    "Gly   {:.2%}\n".format((proteinStrandFrequencyGlobal["Gly"]/proteinStrandFrequencyGlobal["Total"])),
    "Val   {:.2%}\n".format((proteinStrandFrequencyGlobal["Val"]/proteinStrandFrequencyGlobal["Total"])),
    "Ser   {:.2%}\n".format((proteinStrandFrequencyGlobal["Ser"]/proteinStrandFrequencyGlobal["Total"])),
    "\n\n",
    "Total Protiens: "+str(len(translatedData))+"\n\n"
])
results.close()