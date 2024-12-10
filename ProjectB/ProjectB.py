
import math
import inspect

def readFile(filename):
    content=""
    with open(filename,"r") as fin:
        content=fin.read()
    return content

def getTitle(strFile):
    """
    Gets the title from the file which is the first line.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    title : str
        Title of the crystal.
    """
    title = ""
    title = strFile.splitlines()[0] #gets first line of strFile
    return title

def getScalingFactor(strFile):
    """
    Gets the scaling factor which is the 1st element of the 2nd line.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    sfactor : float
        Scaling factor on the 2nd line.
    """
    sfactor = 0.0
    splitlines = strFile.splitlines()[1].split()[0]
    sfactor = float(splitlines) # gets 2nd line and makes the value a float.
    return sfactor

def getElements(strFile):
    """
    Gets the elements which is on the 6th line. 

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    elements : list
        Elements that make up the crystal.
    """
    elements = []
    elements = strFile.splitlines()[5].split() #gets the 6th line and splits
    #it into a list.
    return elements

# It then
#splits that line into just the number of elements if the pocc sites are 
#there too.
def getNumElements(strFile):
    """
    Gets the number of elements which is the 7th line and seperates the
    number of elements from the pocc sites.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    nelements : list
        Number of each element that is in the crystal.
    """
    nelements = []
    split_up  = strFile.splitlines()[6]
    #gets number of elements by iterating over 7th line and seperating pocc.
    for i in range(3):
        nelements.append(int(split_up.split()[i].split('*')[0]))
    return nelements

def getAtomDecorations(strFile):
    """
    Outputs a list of the elements and how may of each there are.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    decor : list
        List of all elements with the number of occurences of each.
    """
    decor = []
    #calls past functions to define variables
    elements  = getElements(strFile)
    nelements = getNumElements(strFile)
    #adds each element to decor by multipying by number of elements
    for i in range(len(elements)): 
        decor = decor + ([elements[i]] * int(nelements[i]))
    return decor

def parseStrFile(strFile):
    """
    returns a list of the title, elements, and the number of elements.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    title : str
        Title of the crystal.
    elements : list
        Elements that make up the crystal.
    nelements : list
        Number of each element that is in the crystal.
    """
    title=""
    elements  = []
    nelements = []
    #calls past functions to define variables that will be returned in a list.
    title = getTitle(strFile)
    elements = getElements(strFile)
    nelements = getNumElements(strFile)
    return title,elements,nelements

#This function replaces one of the elements with a replacement.
def getElementsNew(elements,element,replacement):
    """
    Replaces one element in elements with a replacement element.

    Parameters
    ----------
    elements : list
        Elements that make up the crystal.
    element : str
        Element to be replaced.
    replacement : str
        Replacement element.

    Returns
    -------
    elements_new : list
        Element list with element replaced with replacement.
    """
    elements_new = []
    elements_new = elements
    for i in range(len(elements)):
        if elements_new[i] == element:
            elements_new[i] = replacement
    return elements_new

#This function gets the GCD of the number of elements
def getGCD(lst_int):
    _gcd = math.gcd(*lst_int)
    return _gcd

def getTitleNew(elements_new,nelements_reduced):
    """
    Gets the empirical formula of the molecule.

    Parameters
    ----------
    elements_new : list
        Element list with element replaced with replacement.
    nelements_reduced : list
        Number of each element but reduced by their gcd.

    Returns
    -------
    title_new : str
        New title with the replacement element.
    """
    title_new=""
    #makes new title by iterating through elements new and nelements reduced.
    for i in range(len(elements_new)):  
        title_new += elements_new[i] + str(nelements_reduced[i])
    return title_new

def getReducedComposition(nelements):
    """
    Gets the numbers for the empirical formula using the gcd.

    Parameters
    ----------
    nelements : list
        Number of each element that is in the crystal.

    Returns
    -------
    nelements_reduced : list
        Number of each element that is in the crystal divided by their gcd.
    """
    nelements_reduced = []
    _gcd = getGCD(nelements)
    #iterates through nelements and divides that number by the gcd.
    for i in range(len(nelements)):
        nelements_reduced.append(int((nelements[i]) / _gcd))
    return nelements_reduced

def replaceAtomPositions(lines_new,nelements,elements_new,element,replacement):
    """
    Replaces element with replacement in 
    
    element in par with its replacement. 

    Parameters
    ----------
    lines_new : list
        strFile that has been converted to lines new.
    nelements : list
        not used.
    elements_new : lisst
        not used.
    element : str
        Element that will be replaced.
    replacement : str
        Element that is replacing element.

    Returns
    -------
    lines_new : list
        Lines_new that has had the element replaced with replacement.
    """
    for i in range(8,len(lines_new)):
#iterates through the coordinate lines and replaces element with replacement.
        if lines_new[i].split()[3] == element: 
            lines_new[i] = lines_new[i].replace(element,replacement)
    return lines_new

def replaceElement(strFile,element,replacement):
    """
    Calls other functions in order to replace one element with a 
    replacement element.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.
    element : str
        Element that will be replaced.
    replacement : str
        Element that is replacing element.

    Returns
    -------
    strFile_new : str
        strFile that has had all occurences of element replaced with
        replacement.
    """
    strFile_new = ""
    #initializing variables.
    lines_new = strFile.splitlines()
    title,elements,nelements = parseStrFile(strFile)
    #calling functions to initialize more variables.
    elements_new = getElementsNew(elements,element,replacement)
    nelements_reduced = getReducedComposition(nelements)
    lines_new = replaceAtomPositions(lines_new,nelements,elements_new,
                                     element,replacement)   
    #replacing occurences of element with replacement.
    lines_new[0] = getTitleNew(elements_new,nelements_reduced)
    lines_new[5] = ' '.join(elements)
    strFile_new  = "\n".join(lines_new) 
    return strFile_new

def getPOCCTols(strFile):
    """
    Gets the tolerences for pocc which are on the 2nd line in the 2nd 
    and 3rd positions.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    pocc_tols : list
        LIst of the relative and absolute tolerences.
    """
    pocc_tols = []
    #splits the lines and then gets both the tolerences by iterating in a loop.
    splitlines = strFile.splitlines()[1].split()
    for i in range(2):
        pocc_tols.append(float(splitlines[(i+1)]))
    return pocc_tols

#This function gets the number of elements and the occs which are on the 7th
#line. The # elements and occ are split by '*' and # elements is first 
#and the occs are 2nd in that list. First is int and 2nd is float.
def getNumElementsPOCC(strFile):
    """
    Gets the number of elements and the occs which are on the 7th
    line. The number of elements and occ are split by '*' and the number
    of elements is first and the occs are 2nd in that list.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    nelements : int
        Number of each element in the crystal.
    occs : float
        Partial occupancy values.
    """
    nelements = []
    split_up  = strFile.splitlines()[6]
    #iterates through line 7 and gets the nelements.
    for i in range(3):
        nelements.append(int(split_up.split()[i].split('*')[0]))
    occs = []
    #iterates through line 7 and gets the occs.
    for i in range(3):
        occs.append(float(split_up.split()[i].split('*')[1]))
    return nelements,occs

def coordsAreClose(coord1,coord2,tol=1e-3): 
    """
    Checks if the coordinates are close and returns it to getPOCCSites.

    Parameters
    ----------
    coord1 : list
        Coordinates of the first potential pocc site.
    coord2 : list
        Coordinates of the second potential pocc site.
    tol : float, optional
        Absolute tolerence. The default is 1e-3.

    Returns
    -------
    Boolean
        True or false for if the coords are close or not.

    """
    #iterates through the points in coord1 and coord2 and compares them.
    result = [math.isclose(c1, c2 , abs_tol=tol) for c1, c2 in 
              zip(coord1,coord2)]
    return all(result)

def getPOCCSites(strFile):
    """
    Gets the pairs of POCC sites by seeing the elements that have the
    same coordiantes, uses coords are close.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.

    Returns
    -------
    pocc_sites : list
        Contains the pairs for the elements that are sharing the same
        coordinates.
    """
    #initializing variables
    pocc_sites    = []
    lines_new     = strFile.splitlines()
    iatom_matched = []
    tol           = 1e-3
    nelements     = getNumElements(strFile)
    total         = sum(nelements) + 8

#itereates through all the coordinates and sets coordiante1 as list of floats.
#checks if coord1 has already been checked.
    for i in range(8,total):
        iatom = (i - 8)
        if iatom in iatom_matched:
            break
        coord1 = lines_new[i].split()[0:3]
        coord1floats = [float(y) for y in coord1]
        site_new = True
#itereates through all the coordinates and sets coordiante2 as list of floats.
#Checks if coord2 has already been checked.
        for j in range((i+1),total):
            jatom = (j - 8)
            if jatom not in iatom_matched:
                coord2 = lines_new[j].split()[0:3]
                coord2floats = [float(y) for y in coord2]
#compares coord1 and coord2 and sees if they are close. Adds them to 
#pocc_sites.
                if coordsAreClose(coord1floats,coord2floats,tol):
                    if site_new:
                        iatom_matched.append(iatom)
                        site_new = False
                        pocc_sites.append([iatom,jatom])
                        iatom_matched.append(jatom)   
    return pocc_sites

def getNAtomsPOCC(strFile,supercell_size):
    """
    Returns the number of atoms that there would be if there was no overlap.

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.
    supercell_size : int
        Size of the supercell.

    Returns
    -------
    natoms_ordered : int
        The number of atoms that there would be if there was no overlap.
    """
    natoms_ordered  = 0
    nelements ,occs =  getNumElementsPOCC(strFile)
#iterates through nelements and occs and multiplies them together and adds
#them.
    for (nelements , occs) in zip(nelements ,occs):
        natoms_ordered += float(nelements * occs)
    natoms_ordered *= float(supercell_size)
    return int(natoms_ordered)

def getAtomDecorationsPOCC(strFile,supercell_size):
    """
    Finds the dcorations if the elements were not sharing any 

    Parameters
    ----------
    strFile : str
        String that contains the structure of the crystal.
    supercell_size : int
        Size of the supercell.

    Returns
    -------
    decor : list
        List of all elements with the number of occurences of each.
    """
    #initializing variables
    decor           = []
    total           = []
    elements        = getElements(strFile)
    nelements ,occs = getNumElementsPOCC(strFile)
#Getting summation of nelements and multiplying them by occs. Then multiplying
#it by supercell size.
    for i in range(len(nelements)):
        total.append(float(nelements[i]) * float(occs[i]) * float(supercell_size))
#Adding elements the number of times in total to get decor. 
    for i in range(len(nelements)):
        decor = decor + ([elements[i]] * int(total[i]))
    return decor

def main():
    strFile=readFile("POSCAR_AlNTi2")
    
    ("getTitle()=",getTitle(strFile))
    (getScalingFactor(strFile))
    (getElements(strFile))
    (getNumElements(strFile))
    (getAtomDecorations(strFile))
    (parseStrFile(strFile))
    (getElementsNew(['Al', 'N', 'Ti'],'Ti','C'))
    (getReducedComposition([2, 2, 4]))
    (getTitleNew(['Al', 'C', 'Ti'],[1, 1, 2]))
    (replaceElement(strFile,"N","Na"))
    
    
    
    """
    print("getTitle()=",getTitle(strFile))
    print(getScalingFactor(strFile))
    print(getElements(strFile))
    print(getNumElements(strFile))
    print(getAtomDecorations(strFile))
    print(parseStrFile(strFile))
    print(getElementsNew(['Al', 'N', 'Ti'],'Ti','C'))
    print(getReducedComposition([2, 2, 4]))
    print(getTitleNew(['Al', 'C', 'Ti'],[1, 1, 2]))
    print(replaceElement(strFile,"N","Na"))
    """
    strFile=readFile("PARTCAR_SSeZn")
    """
    print("getTitle()=",getTitle(strFile))
    print(getScalingFactor(strFile))
    print(getPOCCTols(strFile))
    print(getNumElementsPOCC(strFile))
    print(getPOCCSites(strFile))
    print(getNAtomsPOCC(strFile,"4"))
    print(getAtomDecorationsPOCC(strFile,"4"))
    """
if __name__ == "__main__":
    main()
