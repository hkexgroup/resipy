#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 14:23:10 2018
RAM usage
@author: jimmy
"""

#%% function pulled from internet to look at RAM usage
class memoryCheck():
    """Checks memory of a given system"""
 
    def __init__(self):
 
        if os.name == "posix":
            self.value = self.linuxRam()
        elif os.name == "nt":
            self.value = self.windowsRam()
        else:
            print("I only work with Win or Linux :P")
 
    def windowsRam(self):
        """Uses Windows API to check RAM in this OS"""
        kernel32 = ctypes.windll.kernel32
        c_ulong = ctypes.c_ulong
        class MEMORYSTATUS(ctypes.Structure):
            _fields_ = [
                ("dwLength", c_ulong),
                ("dwMemoryLoad", c_ulong),
                ("dwTotalPhys", c_ulong),
                ("dwAvailPhys", c_ulong),
                ("dwTotalPageFile", c_ulong),
                ("dwAvailPageFile", c_ulong),
                ("dwTotalVirtual", c_ulong),
                ("dwAvailVirtual", c_ulong)
            ]
        memoryStatus = MEMORYSTATUS()
        memoryStatus.dwLength = ctypes.sizeof(MEMORYSTATUS)
        kernel32.GlobalMemoryStatus(ctypes.byref(memoryStatus))
 
        return int(memoryStatus.dwTotalPhys/1024**2)
 
    def linuxRam(self):
        """Returns the RAM of a linux system"""
        totalMemory = os.popen("free -m").readlines()[1].split()[1]
        return int(totalMemory)