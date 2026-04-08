import numpy as np

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
        self.primitiveFunctionList = primitiveFunctionList
    
    def evaluate(self, q: float) -> float:
        functionOutput: float = 1.0
        for i in range(self.numberOfFunctions):
            functionOutput *= self.primitiveFunctionList[i](q)
        return functionOutput
        
def readBasicFunctions(basicFunctionInputFile: str, masses : np.ndarray) -> dict:
    with open(basicFunctionInputFile) as f:
        basicFunctionInput: str = f.read()
    basicFunctionsList: dict = {}
    numberOfModes: int = basicFunctionInput.lower().count("mode")
    basicFunctionInputLines: list = basicFunctionInput.split("\nEND")[0].split("\n")
    basicFunctionMainHeader: str = basicFunctionInputLines[0].lower()
    basicFunctionInputLines = basicFunctionInputLines[1:]
    modeHeaderIndex: int = 0
    if "poten" in basicFunctionMainHeader:
        print("die katze ist traurig")
        # for i in range(numberOfModes):
        #     modeFunctionList = {} # New list of functions for mode
        #     # modeFunctionList[0] = basicFunction("0 1 0 I 1 1") # Function 0 is always 1!
        #     numberOfLinesForMode: int = int(basicFunctionInputLines[modeHeaderIndex].split()[-1])
        #     lineBeingRead: int = 1
        #     functionIndexCounter: int = 0
        #     while lineBeingRead <= numberOfLinesForMode:
        #         basicFunctionInputLine: str = basicFunctionInputLines[modeHeaderIndex + lineBeingRead]
        #         basicFunctionInputLineSplit: list = basicFunctionInputLine.split()
        #         functionIndex: int = int(basicFunctionInputLineSplit[0])
        #         if int(basicFunctionInputLineSplit[0]) < 0:
        #             modeFunctionList[functionIndex] = basicFunction(basicFunctionInputLine)
        #         else:
        #             functionLabel: str = basicFunctionInputLineSplit[3].lower()
        #             seriesType: str = "taylor"
        #             orderIndex: int = 2
        #             if seriesType == "taylor":
        #                 orderIndex = 2
        #             else:
        #                 orderIndex = 4
        #             maxOrder: int = int(basicFunctionInputLineSplit[orderIndex])
        #             for j in range(maxOrder+1):
        #                 modeFunctionList[functionIndexCounter] = basicFunction(basicFunctionInputLine, True, j)
        #                 functionIndexCounter += 1
        #         lineBeingRead += 1
        #     basicFunctionsList[i + 1] = modeFunctionList
        #     modeHeaderIndex += lineBeingRead
    else:
        if "mass" in basicFunctionMainHeader:
            basicFunctionsList[0] = 1/masses # Masses are effectively "mode" 0
            for i in range(numberOfModes):
                modeFunctionList = {} # New list of functions for mode
                modeFunctionList[0] = basicFunction("0 1 r 1 1") # Function 0 is always 1!
                numberOfLinesForMode: int = int(basicFunctionInputLines[modeHeaderIndex].split()[-1])
                lineBeingRead: int = 1
                while lineBeingRead <= numberOfLinesForMode:
                    modeFunctionList[lineBeingRead] = basicFunction(basicFunctionInputLines[modeHeaderIndex + lineBeingRead])
                    lineBeingRead += 1
                basicFunctionsList[i] = modeFunctionList
                modeHeaderIndex += lineBeingRead
    return basicFunctionsList

# Labels for kinetic energy operator components
kineticComponents = ["gvib", "grot", "gcor", "pseudo"]

class operatorMap:
    operatorType: str
    operatorRank: int
    componentIndices: dict = {}
    componentCoefficients: dict = {}
    functionIndices: dict = {}
    numberOfModes: int
    containsMass: bool = False

    def __init__(self, operatorMappingInput: str):
        operatorMappingInputLines: list = operatorMappingInput.split("\nEND")[0].split("\n")
        operatorMappingHeader: str = operatorMappingInputLines[0].lower()
        self.operatorType = operatorMappingHeader.split()[1]
        self.containsMass = "mass" in operatorMappingHeader
        self.operatorRank = operatorRankMapping[self.operatorType]
        self.numberOfModes = len(operatorMappingInputLines[2].split()) - self.operatorRank - 3 - self.containsMass
        if self.operatorType == "kinetic":
            pass
        else:
            operatorComponentInputLines: list = operatorMappingInputLines[2 :]
            numberOfComponentTerms: int = len(operatorComponentInputLines)
            self.componentCoefficients[self.operatorType] = np.zeros(numberOfComponentTerms)
            if self.operatorRank > 0:
                self.componentIndices[self.operatorType] = np.zeros((numberOfComponentTerms, self.operatorRank))
            self.functionIndices[self.operatorType] = np.zeros((numberOfComponentTerms, self.numberOfModes + self.containsMass))
            for i in range(numberOfComponentTerms):
                operatorComponentLineSplit: list = operatorComponentInputLines[i].split()
                self.componentCoefficients[self.operatorType][i] = float(operatorComponentLineSplit[self.operatorRank + 2])
                for j in range(self.operatorRank):
                    self.componentIndices[self.operatorType][i, j] = int(operatorComponentLineSplit[j])
                for j in range(self.numberOfModes + self.containsMass):
                    self.functionIndices[self.operatorType][i, j] = int(operatorComponentLineSplit[j + self.operatorRank + 3])
    
    def evaluatePointOfComponent(self, component: str, basicFunctionsList: dict, q):
        coefficients = self.componentCoefficients[component]
        functionIndices = self.functionIndices[component]
        numberOfTerms: int = len(coefficients)
        terms = np.ones(numberOfTerms)
        match self.operatorRank:
            case 0:
                for i in range(numberOfTerms):
                    terms[i] *= coefficients[i]
                    if self.containsMass:
                        terms[i] *= basicFunctionsList[0][functionIndices[i, 0]]
                        for j in range(1, self.containsMass + self.numberOfModes):
                            terms[i] *= basicFunctionsList[j][functionIndices[i, j]].evaluate(q[j])
                    else:
                        for j in range(self.numberOfModes):
                            terms[i] *= basicFunctionsList[j + 1][functionIndices[i, j]].evaluate(q[j])
                return np.sum(terms)
            case 1:
                operatorValue = np.zeros(3)
            case 2:
                match component:
                    case "pseudo":
                        operatorValue = 0.0
                    case "gvib":
                        operatorValue = np.zeros((self.numberOfModes, self.numberOfModes))
                    case "grot":
                        operatorValue = np.zeros((3, 3))
                    case "gcor":
                        operatorValue = np.zeros((self.numberOfModes, 3))