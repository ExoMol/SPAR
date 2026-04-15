import numpy as np
import re

wavenumberConversion = 33.7152537836138732499014581824370 # conversion factor to cm-1

class basicFunction:
    basicFunctionType: str
    numberOfFunctions: int
    primitiveFunctionList: list

    def __init__(self, basicFunctionInputLine: str, isSeries: bool = False, order: int = 0) -> None:
        '''isSeries: Is basic function part of a series expansion
           order: Order for a repeated function type
        '''
        basicFunctionLineSplit: list = basicFunctionInputLine.split()
        primitiveFunctionList: list = []
        
        if isSeries:
            self.numberOfFunctions = int(basicFunctionLineSplit[0])
            readIndex : int = 1
        else:
            self.numberOfFunctions = int(basicFunctionLineSplit[1])
            readIndex : int = 2
        for i in range(self.numberOfFunctions):
            functionString = basicFunctionLineSplit[readIndex].lower()
            match functionString:
                # KEO type primitive functions
                case "r":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: (b*q)**p]
                    readIndex += 3
                case "sin":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: np.sin(b*q)**p]
                    readIndex += 3
                case "cos":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: np.cos(b*q)**p]
                    readIndex += 3
                case "tan":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: np.tan(b*q)**p]
                    readIndex += 3
                case "sec":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: np.cos(b*q)**(-p)]
                    readIndex += 3
                case "csc":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: np.sin(b*q)**(-p)]
                    readIndex += 3
                case "cot":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    primitiveFunctionList += [lambda q, p=p, b=b: np.tan(b*q)**(-p)]
                    readIndex += 3
                # Series type primitive functions
                case "morse":
                    b = float(basicFunctionLineSplit[readIndex + 4])
                    c = float(basicFunctionLineSplit[readIndex + 5])
                    primitiveFunctionList += [lambda q, p=order, b=b, c=c: (1-np.exp(-c*(q-b)))**p]
                    readIndex += 6
                case "power":
                    b = float(basicFunctionLineSplit[readIndex + 4])
                    c = float(basicFunctionLineSplit[readIndex + 5])
                    primitiveFunctionList += [lambda q, p=order, b=b: (q - b)**p]
                    readIndex += 6
                case "cos(q)-cos(q0)":
                    b = float(basicFunctionLineSplit[readIndex + 4])
                    c = float(basicFunctionLineSplit[readIndex + 5])
                    primitiveFunctionList += [lambda q, p=order, b=b: (np.cos(q) - np.cos(b))**p]
                    readIndex += 6
                case "sin(q)-sin(q0)":
                    b = float(basicFunctionLineSplit[readIndex + 4])
                    c = float(basicFunctionLineSplit[readIndex + 5])
                    primitiveFunctionList += [lambda q, p=order, b=b: (np.sin(q) - np.sin(b))**p]
                    readIndex += 6
                case "cos(q0)-cos(q)":
                    b = float(basicFunctionLineSplit[readIndex + 4])
                    c = float(basicFunctionLineSplit[readIndex + 5])
                    primitiveFunctionList += [lambda q, p=order, b=b: (np.cos(b) - np.cos(q))**p]
                    readIndex += 6
                case "sin(q0)-sin(q)":
                    b = float(basicFunctionLineSplit[readIndex + 4])
                    c = float(basicFunctionLineSplit[readIndex + 5])
                    primitiveFunctionList += [lambda q, p=order, b=b: (np.sin(b) - np.sin(q))**p]
                    readIndex += 6
                case "exp":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    c = float(basicFunctionLineSplit[readIndex + 3])
                    primitiveFunctionList += [lambda q, p=p, b=b, c=c: np.exp(c*(q-b)**p)]
                    readIndex += 4
                case "expmin":
                    p = float(basicFunctionLineSplit[readIndex + 1])
                    b = float(basicFunctionLineSplit[readIndex + 2])
                    c = float(basicFunctionLineSplit[readIndex + 3])
                    primitiveFunctionList += [lambda q, p=p, b=b, c=c: np.exp(c*min(-abs(np.pi - q) - (b - np.pi), 0)**p)]
                    readIndex += 4
        self.primitiveFunctionList = primitiveFunctionList
    
    def evaluate(self, q: float) -> float:
        functionOutput: float = 1.0
        for i in range(self.numberOfFunctions):
            functionOutput *= self.primitiveFunctionList[i](q)
        return functionOutput

# Defines series type functions     
seriesTypeFunctions = [
    "morse",
    "power",
    "cos(q)-cos(q0)",
    "sin(q)-sin(q0)",
    "cos(q0)-cos(q)",
    "sin(q0)-sin(q)"
]

def readBasicFunctions(basicFunctionInputFile: str) -> dict:
    with open(basicFunctionInputFile) as f:
        basicFunctionInput: str = f.read()
    basicFunctionsList: dict = {}
    numberOfModes: int = basicFunctionInput.lower().count("mode")
    # print(basicFunctionInput)
    basicFunctionInputLines: list = basicFunctionInput.split("\nEND")[0].split("\n")
    basicFunctionMainHeader: str = basicFunctionInputLines[0].lower()
    basicFunctionInputLines = basicFunctionInputLines[1:]
    modeHeaderIndex: int = 0
    if "poten" in basicFunctionMainHeader or "external" in basicFunctionMainHeader:
        for i in range(numberOfModes):
            modeFunctionList = {} # New list of functions for mode
            modeFunctionList[0] = basicFunction("0 1 r 0 1") # Function 0 is always 1!
            lineNumber = 1
            functionIndex = 1
            try:
                while "mode" not in basicFunctionInputLines[modeHeaderIndex + lineNumber].lower():
                    if basicFunctionInputLines[modeHeaderIndex + lineNumber].split()[1].lower() in seriesTypeFunctions:
                        orderLower = int(basicFunctionInputLines[modeHeaderIndex + lineNumber].split()[2])
                        orderUpper = int(basicFunctionInputLines[modeHeaderIndex + lineNumber].split()[4])
                        for j in range(orderLower, orderUpper + 1):
                            modeFunctionList[functionIndex] = basicFunction(basicFunctionInputLines[modeHeaderIndex + lineNumber], True, j)
                            functionIndex += 1
                        lineNumber += 1
                    else:
                        modeFunctionList[functionIndex] = basicFunction(basicFunctionInputLines[modeHeaderIndex + lineNumber], True)
                        functionIndex += 1
                        lineNumber += 1
            except:
                pass
            basicFunctionsList[i + 1] = modeFunctionList
            modeHeaderIndex += lineNumber
    else:
        for i in range(numberOfModes):
            modeFunctionList = {} # New list of functions for mode
            modeFunctionList[0] = basicFunction("0 1 r 0 1") # Function 0 is always 1!
            numberOfLinesForMode: int = int(basicFunctionInputLines[modeHeaderIndex].split()[-1])
            lineBeingRead: int = 1
            while lineBeingRead <= numberOfLinesForMode:
                modeFunctionList[lineBeingRead] = basicFunction(basicFunctionInputLines[modeHeaderIndex + lineBeingRead])
                lineBeingRead += 1
            basicFunctionsList[i + 1] = modeFunctionList
            modeHeaderIndex += lineBeingRead
    return basicFunctionsList

# Keys for kinetic energy operator dictionaries
kineticComponents = ["gvib", "grot", "gcor", "pseudo"]
kineticComponentsRegex = "gvib|grot|gcor|pseudo" # REGEX expression to split KEO checkpoint file

class kineticMapping:
    kineticCoefficients: dict = {} # Dictionary where each entry is a vector
    kineticComponentIndices: dict = {} # Dictionary where each entry is a matrix of integers (N_L x 2) where N_L is the number of terms for the vibrational/rotational/coriolis/pseudopotential part of the KEO
    kineticBasicFunctionIndices: dict = {} # Dictionary where each entry is a matrix of indices for the basic functions (N_L x [3N-6+1]) with masses and (N_L x [3N-6]) without masses
    massesIncluded: bool
    numberOfModes: int

    def __init__(self, kineticCheckpointFile: str, massesIncluded: bool = True):
        with open(kineticCheckpointFile) as f:
            kineticCheckpointContent: str = f.read().lower()
        kineticCheckpointContentSplit = re.split(kineticComponentsRegex, kineticCheckpointContent)
        kineticVibrationalInput = kineticCheckpointContentSplit[1].split("\n")[1:][:-2]
        kineticRotationalInput = kineticCheckpointContentSplit[2].split("\n")[1:][:-2]
        kineticCoriolisInput = kineticCheckpointContentSplit[3].split("\n")[1:][:-2]
        kineticPseudopotentialInput = kineticCheckpointContentSplit[4].split("\nEnd of kinetic")[0].split("\n")[1:][:-4]
        
        self.massesIncluded = massesIncluded
        self.numberOfModes = len(kineticVibrationalInput[0].split()) - self.massesIncluded - 5
        
        # Parse Vibrational input
        self.kineticComponentIndices["gvib"] = np.zeros((len(kineticVibrationalInput), 2), dtype=int)
        self.kineticCoefficients["gvib"] = np.zeros(len(kineticVibrationalInput))
        self.kineticBasicFunctionIndices["gvib"] = np.zeros((len(kineticVibrationalInput), self.numberOfModes + self.massesIncluded), dtype=int)

        for i in range(len(kineticVibrationalInput)):
            componentLineSplit = kineticVibrationalInput[i].split()
            self.kineticComponentIndices["gvib"][i] = np.array([int(componentLineSplit[0]), int(componentLineSplit[1])])
            self.kineticCoefficients["gvib"][i] = float(componentLineSplit[4])
            for j in range(self.numberOfModes + self.massesIncluded):
                self.kineticBasicFunctionIndices["gvib"][i, j] = int(componentLineSplit[5 + j])

        # Parse Rotational input
        self.kineticComponentIndices["grot"] = np.zeros((len(kineticRotationalInput), 2), dtype=int)
        self.kineticCoefficients["grot"] = np.zeros(len(kineticRotationalInput))
        self.kineticBasicFunctionIndices["grot"] = np.zeros((len(kineticRotationalInput), self.numberOfModes + self.massesIncluded), dtype=int)

        for i in range(len(kineticRotationalInput)):
            componentLineSplit = kineticRotationalInput[i].split()
            self.kineticComponentIndices["grot"][i] = np.array([int(componentLineSplit[0]), int(componentLineSplit[1])])
            self.kineticCoefficients["grot"][i] = float(componentLineSplit[4])
            for j in range(self.numberOfModes + self.massesIncluded):
                self.kineticBasicFunctionIndices["grot"][i, j] = int(componentLineSplit[5 + j])

        # Parse Coriolis input
        self.kineticComponentIndices["gcor"] = np.zeros((len(kineticCoriolisInput), 2), dtype=int)
        self.kineticCoefficients["gcor"] = np.zeros(len(kineticCoriolisInput))
        self.kineticBasicFunctionIndices["gcor"] = np.zeros((len(kineticCoriolisInput), self.numberOfModes + self.massesIncluded), dtype=int)

        for i in range(len(kineticCoriolisInput)):
            componentLineSplit = kineticCoriolisInput[i].split()
            self.kineticComponentIndices["gcor"][i] = np.array([int(componentLineSplit[0]), int(componentLineSplit[1])])
            self.kineticCoefficients["gcor"][i] = float(componentLineSplit[4])
            for j in range(self.numberOfModes + self.massesIncluded):
                self.kineticBasicFunctionIndices["gcor"][i, j] = int(componentLineSplit[5 + j])

        # Parse pseudopotential input
        self.kineticComponentIndices["pseudo"] = np.zeros((len(kineticPseudopotentialInput), 2), dtype=int)
        self.kineticCoefficients["pseudo"] = np.zeros(len(kineticPseudopotentialInput))
        self.kineticBasicFunctionIndices["pseudo"] = np.zeros((len(kineticPseudopotentialInput), self.numberOfModes + self.massesIncluded), dtype=int)

        for i in range(len(kineticPseudopotentialInput)):
            componentLineSplit = kineticPseudopotentialInput[i].split()
            self.kineticCoefficients["pseudo"][i] = float(componentLineSplit[4])
            for j in range(self.numberOfModes + self.massesIncluded):
                self.kineticBasicFunctionIndices["pseudo"][i, j] = int(componentLineSplit[5 + j])

    def evaluate(self, kineticComponentLabel: str, basicFunctions: dict, internalCoordinates: np.ndarray, masses: np.ndarray):
        numberOfComponentCoefficients = len(self.kineticCoefficients[kineticComponentLabel])
        match kineticComponentLabel:
            case "gvib":
                kineticComponent = np.zeros((self.numberOfModes, self.numberOfModes))
            case "grot":
                kineticComponent = np.zeros((3, 3))
            case "gcor":
                kineticComponent = np.zeros((self.numberOfModes, 3))
            case "pseudo":
                kineticComponent = 0.0
        if kineticComponentLabel == "pseudo":
            for i in range(numberOfComponentCoefficients):
                if self.massesIncluded:
                    newTerm = self.kineticCoefficients[kineticComponentLabel][i]/masses[self.kineticBasicFunctionIndices[kineticComponentLabel][i, 0] - 1]
                for j in range(self.numberOfModes):
                    newTerm *= basicFunctions[j + 1][self.kineticBasicFunctionIndices[kineticComponentLabel][i, j + self.massesIncluded]].evaluate(internalCoordinates[j])
                kineticComponent += newTerm 
        else:   
            for i in range(numberOfComponentCoefficients):
                newTerm = self.kineticCoefficients[kineticComponentLabel][i]
                if self.massesIncluded:
                    newTerm /= masses[self.kineticBasicFunctionIndices[kineticComponentLabel][i, 0] - 1]
                for j in range(self.numberOfModes):
                    newTerm *= basicFunctions[j + 1][self.kineticBasicFunctionIndices[kineticComponentLabel][i, j + self.massesIncluded]].evaluate(internalCoordinates[j])
                kineticComponent[self.kineticComponentIndices[kineticComponentLabel][i, 0]-1, self.kineticComponentIndices[kineticComponentLabel][i, 1]-1] += newTerm 
        return kineticComponent*wavenumberConversion
            
class potentialMapping:
    potentialCoefficients: np.ndarray
    potentialBasicFunctionIndices: np.ndarray
    massesIncluded: bool
    numberOfModes: int
    numberOfTerms: int

    def __init__(self, potentialCheckpointFile: str, massesIncluded: bool = False):
        with open(potentialCheckpointFile) as f:
            potentialCheckpointContent: str = f.read().lower()
        potentialInput = potentialCheckpointContent.split("\nEnd of potential")[0].split("\n")[1:][:-3]
        
        self.massesIncluded = massesIncluded
        self.numberOfModes = len(potentialInput[0].split()) - self.massesIncluded - 3
        self.numberOfTerms = len(potentialInput)

        self.potentialCoefficients = np.zeros(self.numberOfTerms)
        self.potentialBasicFunctionIndices = np.zeros((self.numberOfTerms, self.numberOfModes), dtype=int)

        for i in range(self.numberOfTerms):
            potentialLineSplit = potentialInput[i].split()
            self.potentialCoefficients[i] = float(potentialLineSplit[2])
            for j in range(self.numberOfModes):
                self.potentialBasicFunctionIndices[i, j] = int(potentialLineSplit[3 + j])
    
    def evaluate(self, basicFunctions: dict, internalCoordinates: np.ndarray):
        potential: float = 0.0
        for i in range(self.numberOfTerms):
            newTerm = self.potentialCoefficients[i]
            for j in range(self.numberOfModes):
                newTerm *= basicFunctions[j + 1][self.potentialBasicFunctionIndices[i, j] + self.massesIncluded].evaluate(internalCoordinates[j])
            potential += newTerm
        return potential
            