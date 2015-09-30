 #!/usr/bin/env python
 
import os,sys
    
def parseBool(str):
    str = str.lower()
    if str=="true" or str=="yes" or str=="on":
        return True
    else:
        return False
