#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 10:45:49 2020

@author: Sina
"""

import platform
from os.path import expanduser
import os

OS = platform.system()
if OS == 'Windows':
    import winreg
    
#%%
class Settings(object):
    def __init__(self):
        """Populate self.param for specific settings"""
        self.param = {} # Dict of settings

    
    def setSetting(self):
        """creates the setting string
        
        Returns
        -------
        string of settings value to be stored locally
        """
        if self.param != {}:
            return '\n'.join('{} = {}'.format(str(key), str(value)) for key, value in self.param.items())
    
    
    def readSetting(self, settings):
        """Reads the retrieved local settings string and populates self.param
        
        Parameters
        ----------
        settings : str
            Raw settings string
        """
        settingList = settings.split('\n')
        for val in settingList:
            key = val.split(' = ')[0]
            value = val.split(' = ')[-1]
            self.param[key] = value
            
    
    # for windows only
    def setReg(self, name, value, REG_PATH=r'SOFTWARE\ResIPy'):
        """Sets registry (key and value), works fine even on multiple user computers. No need to change HKEY_LOCAL_MACHINE
        
        Parameters
        ----------
        name : str
            registry key name
        value : str
            registry key value
        REG_PATH : str, optional
            registry key path in Windows registry
        """
        try:
            winreg.CreateKeyEx(winreg.HKEY_CURRENT_USER, REG_PATH, 0)
            registry_key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, REG_PATH, 0, 
                                           winreg.KEY_WRITE)
            winreg.SetValueEx(registry_key, name, 0, winreg.REG_SZ, value)
            winreg.CloseKey(registry_key)
            
        except WindowsError:
            return None
    
    # for windows only
    def getReg(self, name, REG_PATH=r'SOFTWARE\ResIPy'):
        """Retrieves registry value
        
        Parameters
        ----------
        name : str
            registry key name
        value : str
            registry key value
        REG_PATH : str, optional
            registry key path in Windows registry
        
        Returns
        -------
        value : str
            windows registry key value
        """
        try:
            registry_key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, REG_PATH, 0)
            value, regtype = winreg.QueryValueEx(registry_key, name)
            winreg.CloseKey(registry_key)
            return value
        
        except WindowsError:
            return None   
        
    # for windows only    
    def delReg(self, name, REG_PATH=r'SOFTWARE\ResIPy'):
        """Deletes registry value in a key - not the key
        
        Parameters
        ----------
        name : str
            registry key name
        value : str
            registry key value
        REG_PATH : str, optional
            registry key path in Windows registry
        """
        try:
            registry_key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, REG_PATH, 0, 
                                           winreg.KEY_WRITE)
            winreg.DeleteValue(registry_key, name)
            winreg.CloseKey(registry_key)

        except WindowsError:
            return None 
         
        
    def genLocalSetting(self):
        """saves settings on windows, MacOS and Linux"""
        
        try:
            settingRaw = self.setSetting()
            if OS == 'Windows':
                self.setReg('settings', settingRaw)
                
            elif OS == 'Darwin':
                userDir = expanduser('~/Library/Preferences/')
                with open(os.path.join(userDir,'com.resipy.plist'), 'w') as settings_file:
                    settings_file.write(settingRaw)
            
            elif OS == 'Linux':
                userDir = expanduser('~/.local/share/')
                with open(os.path.join(userDir,'resipy.stt'), 'w') as settings_file:
                    settings_file.write(settingRaw)
            print('Local settings saved!')
            
        except:
            print('ERROR!') # expect the unexpected!
            print('!! SETTINGS NOT SAVED !!')
            pass

    
    def retLocalSetting(self):
        """retrieves local settings file and prepares settings"""
        
        try:
            if OS == 'Windows':
                settingRaw = self.getReg('settings')
                
            elif OS == 'Darwin':
                userDir = expanduser('~/Library/Preferences/')
                with open(os.path.join(userDir,'com.resipy.plist'), 'r') as settings_file:
                    settingRaw = settings_file.read()
            
            elif OS == 'Linux':
                userDir = expanduser('~/.local/share/')
                with open(os.path.join(userDir,'resipy.stt'), 'r') as settings_file:
                    settingRaw = settings_file.read()
            
            self.readSetting(settingRaw)
            print('Local settings retrieved!')
            return True
            
        except:
            print('Local settings not found!')
            return None