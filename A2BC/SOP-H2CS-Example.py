#### Legacy Example
import numpy as np

toRadians = np.pi/180
wavenumberConversion = 33.7152537836138732499014581824370 # conversion factor to cm-1
# Case of H2CS
masses = np.array([12.00000000, 31.9720707, 1.00782503223, 1.00782503223])
qVector = np.array([1.608952, 1.086848, 1.086848, 121.750*toRadians, 121.750*toRadians, np.pi])
numberOfModes = len(qVector)

# Read inputs
basicFunctionsInputFile = "combination.chk"
with open(basicFunctionsInputFile) as f:
    basicFunctionsInput = f.read()

kineticCheckpointFile = "kinetic.chk"
with open(kineticCheckpointFile) as f:
    kineticInput = f.read()

basicFunctionsInputLines = basicFunctionsInput.split("\nEND")[0].split("\n")
basicFunctionsInputLines = basicFunctionsInputLines[1:]

elementaryFunctionMap = {
    "R": lambda x : x,
    "COS": np.cos,
    "SIN": np.sin,
    "TAN": np.tan,
    "SEC": lambda x : 1/np.cos(x),
    "CSC": lambda x : 1/np.sin(x),
    "COT": lambda x : 1/np.tan(x)
}

def basicFunctionTemplate(basicFunctionLine: str):
    '''Template for the basic functions, input is a string as defined
    in the basic functions block example of the paper. Output is a 
    function f_p_k.'''
    basicFunctionLineSplit = basicFunctionLine.split()
    if basicFunctionLineSplit[1] == "2":
        f1 = elementaryFunctionMap[basicFunctionLineSplit[2].upper()]
        n1 = float(basicFunctionLineSplit[3])
        a1 = float(basicFunctionLineSplit[4])
        f2 = elementaryFunctionMap[basicFunctionLineSplit[5].upper()]
        n2 = float(basicFunctionLineSplit[6])
        a2 = float(basicFunctionLineSplit[7])
        return lambda q : f1(a1*q)**n1*f2(a2*q)**n2
    else:
        f1 = elementaryFunctionMap[basicFunctionLineSplit[2].upper()]
        n1 = float(basicFunctionLineSplit[3])
        a1 = float(basicFunctionLineSplit[4])
        return lambda q : f1(a1*q)**n1
    
basicFunctionsList = {}
basicFunctionsList[0] = 1/masses # Reciprocal masses act as mode 0
modeBeingRead = 1
headerIndex = 0
while modeBeingRead <= numberOfModes:
    modeFunctionList = {} # New list of functions for mode
    modeFunctionList[0] = lambda q : 1 # Function 0 is always 1!
    numberOfFunctionsForMode = int(basicFunctionsInputLines[headerIndex].split()[-1])
    functionIndex = 1
    while functionIndex <= numberOfFunctionsForMode:
        modeFunctionList[functionIndex] = basicFunctionTemplate(basicFunctionsInputLines[headerIndex + functionIndex])
        functionIndex += 1
    basicFunctionsList[modeBeingRead] = modeFunctionList
    headerIndex += functionIndex
    modeBeingRead += 1

# Divide kinetic input blocks
keywords = ["Gvib\n", "Grot\n", "Gcor\n", "pseudo\n"]
padding = "\n987654321    0    0    0    0.00000000E+00   0"
for i in range(numberOfModes):
    padding += "    0"
padding += "\n"
kineticInputSplit = kineticInput.split(padding)[:4]
kineticVibrationalInput = kineticInputSplit[0].split(keywords[0])[1]
kineticRotationalInput = kineticInputSplit[1].split(keywords[1])[1]
kineticCoriolisInput = kineticInputSplit[2].split(keywords[2])[1]
kineticPseudoPotentialInput = kineticInputSplit[3].split(keywords[3])[1]

GMatrixVibrational = np.zeros((numberOfModes, numberOfModes))
GMatrixRotational = np.zeros((3, 3))
GMatrixCoriolis = np.zeros((numberOfModes, 3))
pseudoPotential = 0.0

kineticPseudoPotentialInputLines = kineticPseudoPotentialInput.split("\n")
for line in kineticPseudoPotentialInputLines:
    kineticLineSplit = line.split()
    newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
    for k in range(1, numberOfModes+1):
        newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
    pseudoPotential += newTerm


kineticVibrationalInputLines = kineticVibrationalInput.split("\n")
for line in kineticVibrationalInputLines:
    kineticLineSplit = line.split()
    newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
    for k in range(1, numberOfModes+1):
        newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
    GMatrixVibrational[int(kineticLineSplit[0])-1, int(kineticLineSplit[1])-1] += newTerm

kineticRotationalInputLines = kineticRotationalInput.split("\n")
for line in kineticRotationalInputLines:
    kineticLineSplit = line.split()
    newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
    for k in range(1, numberOfModes+1):
        newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
    GMatrixRotational[int(kineticLineSplit[0])-1, int(kineticLineSplit[1])-1] += newTerm

kineticCoriolisInputLines = kineticCoriolisInput.split("\n")
for line in kineticCoriolisInputLines:
    kineticLineSplit = line.split()
    newTerm = float(kineticLineSplit[4])*basicFunctionsList[0][int(kineticLineSplit[5])-1]
    for k in range(1, numberOfModes+1):
        newTerm *= basicFunctionsList[k][int(kineticLineSplit[5 + k])](qVector[k-1])
    GMatrixCoriolis[int(kineticLineSplit[0])-1, int(kineticLineSplit[1])-1] += newTerm

pseudoPotential *= wavenumberConversion
GMatrixVibrational *= wavenumberConversion
GMatrixRotational *= wavenumberConversion
GMatrixCoriolis *= wavenumberConversion
print(pseudoPotential)
print(GMatrixVibrational)
print(GMatrixRotational)
print(GMatrixCoriolis)