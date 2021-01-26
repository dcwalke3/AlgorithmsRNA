'''
By: Dakota Walker
Class: Algorithms
Assignment: Protein Translation
'''

#Gets file containing the Raw DNA code.
with open('C:\\Users\dakot\Documents\Algorithms\RNA-to-DNA\RawProtein.txt', 'r') as file:
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

print(protienCountGlobal)


# Seperates things into a group of 3.
n=3
seperatedData = [newData[i:i+n] for i in range(0, len(newData), n)]



#Used to count protein Frequency per protein.
proteinFrequencyCounter = {"A":0, "C":0, "G":0, "U":0, "Total":0}

#Counts the number of certain types of protein
protienStrandFrequencyGlobal = {"Asp":0, "Asn":0, "Arg":0, "Tyr":0, "Gln":0,
                                "Trp":0, "Glu":0, "His":0, "Cys":0, "Pro":0,
                                "Ile":0, "STOP":0, "Thr":0, "Met":0, "Leu":0,
                                "Ala":0, "Gly":0, "Val":0, "Ser":0, "Phe":0,
                                "Lys":0 ,"Total":0}

#Each number should equal a protien Frequency Dictionary
NumbersCounted = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0}


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
                        "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAC":"Glu",
                        "GGU":"Gly", "GGC":"Gly", "GGG":"Gly", "GGA":"Gly",}

translatedData = []

def Conversion():
    for item in seperatedData:
        for key in ProteinConversionKey:
            if(item == key):
                translatedData.append(ProteinConversionKey[key]+"-")


Conversion()
translatedDataString = " "
translatedDataString = translatedDataString.join(translatedData)
translatedDataString = translatedDataString.replace(" ", "").replace("STOP-", "STOP\n")

CountList = translatedDataString.split("\n")
print(CountList)
for c in range(len(CountList)):
    CountList[c] = CountList[c].split("-")

for c in range(len(CountList)):
    for z in range(len(CountList[c])):
        protienStrandFrequencyGlobal[CountList[c][z]]+=1
        protienStrandFrequencyGlobal["Total"]+=1
print(CountList)
print (protienStrandFrequencyGlobal)
