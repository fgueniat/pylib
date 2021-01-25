import fileinput
import re


def Update_field_Control(fname,field,value):
    '''
        update a field in the controlDict in open foam
        fname is the path to the controlDict
        field is the name of the field (e.g. 'startTime')
        value is the associated value
    '''
    for line in fileinput.input(fname,inplace=True):
        if field in line:
            # remove extra blanks:
            line = re.sub(' +',' ',line)
            # make a list for manipulation
            line_split = line.split(' ')
            if line_split[0]==field:# the line has to start with the field !!
                #the first item is the value, so let's replace it
                line_split[1] = str(value) + ';\n'
                line = ' '.join(line_split)
        #the file is entirely rewritten. new lines are in place
        print(line, end='')
    fileinput.close() 


    
def Update_field_UBC(fname,field,values):
    '''
        update the BC field in the boundaryField in open foam
        fname is the path to the field U
        field is the name of the field (e.g. 'movingWall')
        value is a dictionarry of the the associated subfield and values.
    '''
    start_write=False

    for line in fileinput.input(fname,inplace=True):
        if field in line:
            start_write=True
        if '}' in line:
            start_write = False
        if start_write is True:
            for value in values.keys():
                if value in line:
                    # remove extra blanks:
                    line = re.sub(' +',' ',line)
                    # make a list for manipulation
                    line_split = line.split(' ')
                    line_split  = [value,values[value]+';\n']
                    line =  '        ' + '   '.join(line_split) 
        print(line, end='')
    fileinput.close()

