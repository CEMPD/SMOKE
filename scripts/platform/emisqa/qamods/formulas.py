from __future__ import print_function
from builtins import range
from re import match
from qamods.species_array import SpeciesArray

def calc_form(outDict, formK, formNK, ignoreSpec, verbosity):
    '''
    Calculates output species from formulas
    '''
    rmSpec = []
    
    for x in range(2):
        if x == 0:
            formList = formK
            keep = True
        else:
            formList = formNK
            keep = False

        if not formList:
            continue

        formList = [formula.strip() for formula in formList.split(',')]
    
        # Loop over the formulas    
        for eq in formList:

            # Get the output species name
            outSpec = eq.split('=')[0]
            if outSpec in list(outDict.keys()):
                raise ValueError('Output species name %s already exists in species list')
            formula = eq.split('=')[1]

            # Parse the formula
            pol = ''
            formOut = ''
            for y,c in enumerate(formula):
                if match('[A-Z]', c):
                    pol += c
                elif match('[0-9\_\-]', c) and pol !='':
                    # Put numbers in pol names, keeping in mind that pols can't start with a number
                    pol += c
                else:
                    if pol:
                        if pol not in list(outDict.keys()):
                            if ignoreSpec:
                                formOut += '0'
                                print('Warning: Input species %s does not exist.  Replacing with 0 in formula for %s.' %(pol, outSpec))
                            else:
                                raise ValueError('Cannot calculate %s.  Input species %s does not exist.\nMake sure that species is specified after -s or use -a.' %(outSpec,pol))
                        else:
                            formOut += 'outDict[\'%s\']()' %pol
                            if not keep and pol not in rmSpec:
                                rmSpec.append(pol)
                        pol = ''
                    formOut += c

                # Append the last pollutant of the formula if there is one
                if y == (len(formula) - 1) and pol:
                    if pol not in list(outDict.keys()):
                        if ignoreSpec:
                            formOut += '0'
                            print('Warning: Input species %s does not exist.  Replacing with 0 in formula for %s.' %(pol, outSpec))
                        else:
                            raise ValueError('Cannot calculate %s.  Input species %s does not exist.\nMake sure that species is specified after -s or use -a.' %(outSpec,pol))
                    else:
                        formOut += 'outDict[\'%s\']()' %pol
                        if not keep and pol not in rmSpec:
                            rmSpec.append(pol)
 
            if verbosity:
                print('Calculating %s = %s' %(outSpec, formOut))

            # Calculate the formula if it at least one input is present, otherwise don't write
            if 'outDict' in formOut:
                outDict[outSpec] = SpeciesArray(eval(formOut), outSpec)

    # Delete output pollutants marked for removal
    for pol in rmSpec:
        del outDict[pol]

    return outDict


