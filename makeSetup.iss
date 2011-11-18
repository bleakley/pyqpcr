; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{AB0FCAB4-3161-4336-99F1-77E465DE2A79}
AppName=pyQPCR
AppVersion=0.9
AppVerName=pyQPCR - 0.9
AppContact=tgastine[at]users.sourceforge.net
AppPublisher=pyqpcr.sourceforge.net
AppPublisherURL=http://pyqpcr.sourceforge.net
AppSupportURL=https://sourceforge.net/projects/pyqpcr/forums
AppUpdatesURL=http://pyqpcr.sourceforge.net
AppComments=pyQPCR is an open-source software. It can be used to compute qPCR analysis.
DefaultDirName={pf}\pyQPCR
DefaultGroupName=pyQPCR
LicenseFile=C:\Documents and Settings\Thomas\Bureau\pyqpcr\COPYING.txt
OutputBaseFilename=setup
Compression=lzma
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "french"; MessagesFile: "compiler:Languages\French.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "C:\Documents and Settings\Thomas\Bureau\pyqpcr\dist\qpcr.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Documents and Settings\Thomas\Bureau\pyqpcr\dist\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\pyQPCR"; Filename: "{app}\qpcr.exe"
Name: "{commondesktop}\pyQPCR"; Filename: "{app}\qpcr.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\qpcr.exe"; Description: "{cm:LaunchProgram,pyQPCR}"; Flags: nowait postinstall skipifsilent

