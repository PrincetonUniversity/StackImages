#!/usr/bin/env python3

from numpy import *
from astropy.utils.data import download_file
from astropy.io import fits
import wx  
import platform

global folder
folder = '.'

# cached downloaded files are here...
#     ~/.astropy/cache/download/url

def GetFit(title): 
   global tc, tc1, tc2, panel, duh
   
   frame = wx.Frame(None, title=title, size=(800,600))
   frame.SetSize((500,700))
   
   panel = wx.Panel(frame)
   sizer = wx.GridBagSizer(0,0)
           
   text = wx.StaticText(panel, label = "Website:") 
   sizer.Add(text, pos = (0, 0), flag = wx.ALL, border = 5)
           
   tc = wx.TextCtrl(panel) 
   sizer.Add(tc, pos = (0, 1), span = (1, 3), flag = wx.EXPAND|wx.ALL, border = 5) 
    
   text1 = wx.StaticText(panel, label = "File list:") 
   sizer.Add(text1, pos = (1, 0), flag = wx.ALL, border = 5) 
           
   tc1 = wx.TextCtrl(panel,style = wx.TE_MULTILINE) 
   sizer.Add(tc1, pos = (1,1), span = (1, 3), flag = wx.EXPAND|wx.ALL, border = 5) 
   sizer.AddGrowableRow(1) 

   buttonDownload = wx.Button(panel, label = "Download to:") 
   panel.Bind(wx.EVT_BUTTON, pick_folder, buttonDownload)
   sizer.Add(buttonDownload, pos = (2, 0),flag = wx.ALL, border = 5) 
           
   tc2 = wx.StaticText(panel) 
   sizer.Add(tc2, pos = (2,1), span = (1, 3), flag = wx.EXPAND|wx.ALL, border = 5) 
    
   buttonDemo = wx.Button(panel, label = "Click here for an example") 
   panel.Bind(wx.EVT_BUTTON, do_demo, buttonDemo)
   sizer.Add(buttonDemo, pos = (4, 0), span = (1,2), flag = wx.ALL, border = 5) 
   
   buttonOk = wx.Button(panel, label = "Download Files") 
   panel.Bind(wx.EVT_BUTTON, get_them, buttonOk)
   sizer.Add(buttonOk, pos = (4, 2),flag = wx.ALL, border = 5) 
   
   spaces = "                                                                  "
   duh = wx.StaticText(panel, label = spaces)
   sizer.Add(duh, pos = (5, 1), flag = wx.ALL, border = 5)
           
   panel.SetSizerAndFit(sizer)
   
   frame.Centre() 
   frame.Show()      
         
def pick_folder(evt):
   global tc, tc1, tc2, panel, folder, duh

   dlg = wx.DirDialog(panel,"Choose a directory", style=wx.DD_NEW_DIR_BUTTON)
   if dlg.ShowModal() == wx.ID_OK:
       folder = dlg.GetPath()
       tc2.SetLabel(folder+' ')
           
def get_them(evt):
   global tc, tc1, tc2, folder

   url = tc.GetValue()
   files = tc1.GetValue()

   files_list = [y for y in (x.strip() for x in files.splitlines()) if y]

   N = len(files_list)

   duh.SetLabel("Going after the Fits files")

   image_list = [ download_file(url+'/'+file, cache=False, show_progress=False ) for file in files_list ]
   print('OS = '+platform.system())
   if (platform.system() == 'Windows'): slash = '\\'
   else:                                slash = '/'
   for k in range(N):
       duh.SetLabel('Getting '+files_list[k])
       image = fits.getdata(image_list[k])
       image = array(image)
       hdul = fits.open(image_list[k])   # FYI... HDUL = Header Data Unit List
       duh.SetLabel('Saving '+files_list[k])
       fits.writeto(folder+slash+files_list[k], image, hdul[0].header, overwrite=True)

def do_demo(evt):
   global tc, tc1, tc2, panel, folder

   tc.SetValue('https://vanderbei.princeton.edu/fits_files/m27/')
   tc1.SetValue("""m27-0001Ha.fit
m27-0002Ha.fit
m27-0003Ha.fit
m27-0004Ha.fit
m27-0005Ha.fit
m27-0006Ha.fit
m27-0007Ha.fit
m27-0008Ha.fit
m27-0001OIII.fit
m27-0002OIII.fit
m27-0003OIII.fit
m27-0004OIII.fit
m27-0005OIII.fit
m27-0006OIII.fit
m27-0007OIII.fit
m27-0008OIII.fit""")
           
app = wx.App() 
GetFit('Download Astro Fits Files') 
app.MainLoop()
