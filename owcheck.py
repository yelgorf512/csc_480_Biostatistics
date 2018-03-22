#owcheck - module with function to check with user before overwriting pre-existing file
import os

def overwriteFile(filename) :
    if os.path.isfile(filename) :
        yn_input = input(filename + " already exists. Do you want to overwrite it? (Y/N) ")
        if yn_input.upper() == "Y" :
            print("Overwriting...")
            return True
        else :
            print("Okay, exiting without writing...")
            exit(-1)
    else :
        return True