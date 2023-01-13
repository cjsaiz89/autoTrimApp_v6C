# autoTrimApp_v6C
SeaGlider parameter auto trimmer - FTP Online version

Summary of use:

1- Select Working directory in CONFIG tab.
2- Select glider, No. of dives, and last .nc to process. Hit Start.
3- Go to the SFTP tab and connect to the server. Select the glider you're working on.
4- Go to the CMDFILE tab and import the remote cmdfile, then hit Create to display the suggested cmdfile.
5- Copy suggested cmdfile to the left. Edit whatever else you need in the cmdfile and hit Export and accept. The file will be uploaded with the name given in the CONFIG tab.
5- Go to the RESULTS tab and repeat with the next glider. 

1_CONFIG TAB
![4-CONFIG](https://user-images.githubusercontent.com/89260258/212207492-6447d5da-bd56-4290-9bb7-5420002cf61c.PNG)

Once the app starts, go to the CONFIG tab and select your working directory.
This directory will contain the sgGGG folders with its own NCDF files inside of each if you have any. If not, when you download anything it'll automatically create such folders.
Set the parameters to get the suggested calculated values, and change the filter to see only the folders you're interested in from the server.
Save the configuration so next 
time you start the app it'll auto load.

2_RESULTS TAB
![1-Results](https://user-images.githubusercontent.com/89260258/212207500-cc4f0184-d5b7-49bf-9040-f7cb93cd06b4.PNG)


If you have sgGGG folders already in the local working directory, they'll be listed here.
Select the folder and .nc file from which you want to process the dives backwards. Ex.: choosing p6100666.nc and No. of dives 3, will process and average dives 666, 665, and 664.
Press start to see the results. For the SMCC, make sure at least for one of the selected dives is not a NaN.
You can update and automatically download the latest x dives if you're already connected to the server. (x = No. of dives)

3_SFTP TAB
![2-SFTP](https://user-images.githubusercontent.com/89260258/212207518-242d2802-1987-4985-b628-4ab392d4622b.PNG)


Enter user and password to connect and press "Remember me" to save your credentials for next time you use the app.
Here you can download the files you want and the filter in the config tab applies to this dropdown menu.
Press Select to set your current directory to the glider you want to work with. This will be used in the CMDFILE tab.

4_CMDFILE TAB
![3-CMDFILE](https://user-images.githubusercontent.com/89260258/212207588-2e076d06-de68-4f8b-b4d0-5166a208d502.PNG)

On the left you have a whiteboard with what's downloaded and what's uploaded. Whatever you write here will be saved with the filename defined in the CONFIG tab, and sent when you press "Export" to the remote folder selected in the sftp tab.
On the right, you have the suggested cmdfile with the implied changes highlighted, including the config you used. In parentheses you'll see the difference with the parameters that were used in the last processed dive (not the current cmdfile on the basestation).
The imported remote cmdfile file is used as a base to create the implied cmdfile. It'll copy all the parameters in that file, but will only change those that are configured to change (NDIVES, STOPT, SMCC, CVBD, PITCH and ROLL).
If you're happy with the recommended values, press "Copy" to copy the files to the left box. You can still adjust the textbox on the left since it's the final version that will be uploaded when you hit Export.
Remember to change the filename in the config file from cmdfile_upload to cmdfile if you really want to replace the cmdfile
You need to have the cmdfile imported and the dives processed before being able to Create the suggested cmdfile. Also, both sgGGG labels in yellow need to match in order to create the file, copy or export, otherwise it won't let you do it. 
At the bottom, you can see the last time the remote cmdfile imported was modified in the server.
