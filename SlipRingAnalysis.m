% SLIP RING ANALYSIS GUI v.1.0
% This Code Is Written By Fayyaz Nosouhi
% University of Tehran, 15, August 2018

function varargout = SlipRingAnalysis(varargin)
% SLIPRINGANALYSIS MATLAB code for SlipRingAnalysis.fig
%      SLIPRINGANALYSIS, by itself, creates a new SLIPRINGANALYSIS or raises the existing
%      singleton*.
%
%      H = SLIPRINGANALYSIS returns the handle to a new SLIPRINGANALYSIS or the handle to
%      the existing singleton*.
%
%      SLIPRINGANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLIPRINGANALYSIS.M with the given input arguments.
%
%      SLIPRINGANALYSIS('Property','Value',...) creates a new SLIPRINGANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SlipRingAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SlipRingAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SlipRingAnalysis

% Last Modified by GUIDE v2.5 06-Dec-2018 01:24:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SlipRingAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @SlipRingAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SlipRingAnalysis is made visible.
function SlipRingAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)

%Material properties---------------------------
handles.E_SLR=str2double('');
handles.v_SLR=str2double('');
handles.p_SLR=str2double('');
handles.alf_SLR=str2double('');
handles.K_SLR=str2double('');
handles.n_plastic_SLR=str2double('');
handles.table_SLR=[];

handles.E_SHFT=str2double('');
handles.v_SHFT=str2double('');
handles.p_SHFT=str2double('');
handles.alf_SHFT=str2double('');
handles.K_SHFT=str2double('');
handles.n_plastic_SHFT=str2double('');
handles.table_SHFT=[];

handles.E_INS=str2double('');
handles.v_INS=str2double('');
handles.p_INS=str2double('');
handles.alf_INS=str2double('');
handles.K_INS=str2double('');


set(handles.e_slr,'string','');
set(handles.v_slr,'string','');
set(handles.p_slr,'string','');
set(handles.k_slr,'string','');
set(handles.alf_slr,'string','');
set(handles.n_plastic_slr,'string','');

set(handles.e_shft,'string','');
set(handles.v_shft,'string','');
set(handles.p_shft,'string','');
set(handles.k_shft,'string','');
set(handles.alf_shft,'string','');
set(handles.n_plastic_shft,'string','');

set(handles.e_ins,'string','');
set(handles.v_ins,'string','');
set(handles.p_ins,'string','');
set(handles.k_ins,'string','');
set(handles.alf_ins,'string','');

set(handles.table_slr, 'data', {'', ''});
set(handles.table_shft, 'data', {'', ''});
handles.save_mat_index=0;

%Dimension------------------------------------
handles.Dimen1=str2double('');
handles.Dimen2=str2double('');
a='[mm]';
b='[-]';
c={'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','';
    a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, b, b};
set(handles.dimen,'data', c');
handles.save_dim_index=0;
handles.n_old=-1;
%temperature----------------------------------

set(handles. t_teeth, 'ColumnName', 'T[C]');

set(handles. t_teeth, 'ColumnEditable', true);
set(handles. t_teeth, 'data', {'','','','';'','','','' }');

set(handles. theta, 'data', {'','','','','','','','','';'','','','','','','','','' }');
set(handles. theta, 'ColumnEditable', true);
handles.save_temp_index=0;
handles.table_temp1=str2double({'','','','','','','','','';'','','','','','','','',''}');
handles.table_temp2=str2double({'','','','','','','','',''}');
handles.table_teeth2=str2double('');
handles.table_teeth1=str2double('');

% other data--------------------------------
handles.save_other_index=0;
handles.table_therm_press=str2double('');
handles.therm_press_PNT=str2double('');
handles.fric_COF=str2double('');
handles.bolt_LOAD=str2double('');
handles.mass_CONEC=str2double('');
handles.rot_SPEED=str2double('');
handles.mesh_INS=str2double('');
handles.mesh_SLR=str2double('');
handles.mesh_SHFT=str2double('');
handles.I11=str2double('');
handles.I22=str2double('');
handles.I33=str2double('');

set(handles.fric_cof,'string','');
set(handles.mass_conec,'string','');
set(handles.bolt_pre_load,'string','');
set(handles.rot_speed,'string','');
set(handles.mesh_slr,'string','');
set(handles.mesh_ins,'string','');
set(handles.mesh_shft,'string','');
set(handles.thermal_press_point,'string','');
set(handles.thermal_pressure_table, 'data', {'','''',''});
set(handles.i11,'string','');
set(handles.i22,'string','');
set(handles.i33,'string','');
%SN Data--------------------------------------
set(handles.s_n_data,'RowName', {'SlipRing','Shaft','Insulation',});
set(handles.s_n_data,'data', {'','','','','','','','';'','','','','','','','';'','','','','','','','';});
handles.S_N_Data=str2double({'','','','','','','','';'','','','','','','','';'','','','','','','','';});

handles.su_SLR=str2double('');
handles.su_SHFT=str2double('');
handles.su_INS=str2double('');

set(handles.su_slr,'string','');
set(handles.su_shft,'string','');
set(handles.su_ins,'string','');

handles.save_sn_index=0;

%Run------------------------------------------

handles.low_cycle_Fatigue=0;
handles.high_cycle_Fatigue=0;
handles.mesh_Conv=0;

handles.sep_Check=0;

set(handles.low_cycle_fatigue, 'value',0);
set(handles.high_cycle_fatigue, 'value',0);
set(handles.mesh_conv, 'value',0);
set(handles.sep_check, 'value',0);

handles.path='';
handles.model_Name='';
 handles.model_Name='';
handles.save_Dir='';
set(handles.speed_tol,'Enable', 'off')
set(handles.mesh_ref_fac,'Enable', 'off')
set(handles.over_speed,'Enable', 'off')
set(handles.high_cycle_no,'Enable', 'off')
set(handles.low_cycle_no,'Enable', 'off')

handles.high_cycle_No=str2double('');
handles.low_cycle_No=str2double('');
handles.mesh_ref_Fac=str2double('');
handles.over_Speed=str2double('');
handles.speed_Tol=str2double('');

set(handles.high_cycle_no,'string','');
set(handles.low_cycle_no,'string','');
set(handles.over_speed,'string','');
set(handles.mesh_ref_fac,'string','');
set(handles.speed_tol,'string','');
 handles.numCPU=2;
set(handles.num_cpu,'string','2'); 
handles.singular='Linear';


set(handles.model_name,'string','');
set(handles.save_dir,'string','');


set(handles.lin1,'string','11');
set(handles.lin2,'string','16');

set(handles.q1,'string','14');
set(handles.q2,'string','19');
set(handles.q3,'string','24');

set(handles.uipanel_results,'visible', 'off')

handles.LIN1=11;
handles.LIN2=16;
handles.Q1=14;
handles.Q2=19;
handles.Q3=24;
  

set(handles.listbox5, 'string', '');
set(handles.listbox6, 'string', '');

handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SlipRingAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SlipRingAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Material properties
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SLIP RING MAterials PROPERTIES-----------------------------------------------------
function n_plastic_slr_CreateFcn(hObject, eventdata, handles)
function e_slr_CreateFcn(hObject, eventdata, handles)
function v_slr_CreateFcn(hObject, eventdata, handles)
function p_slr_CreateFcn(hObject, eventdata, handles)
function alf_slr_CreateFcn(hObject, eventdata, handles)
function k_slr_CreateFcn(hObject, eventdata, handles)
%
function e_slr_Callback(hObject, eventdata, handles)
   handles.E_SLR=str2double(get(hObject,'String'));
    if isnan(handles.E_SLR) || handles.E_SLR<=0
       errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)

%-------------------------------------------------------------------------
function v_slr_Callback(hObject, eventdata, handles)
    handles.v_SLR=str2double(get(hObject,'String'));
    if isnan(handles.v_SLR) || handles.v_SLR<=0
      errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)

%--------------------------------------------------------------------------
function p_slr_Callback(hObject, eventdata, handles)
    handles.p_SLR=str2double(get(hObject,'String'));
    if isnan( handles.p_SLR) ||  handles.p_SLR<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)


function alf_slr_Callback(hObject, eventdata, handles)
   handles.alf_SLR=str2double(get(hObject,'String'));
    if isnan( handles.alf_SLR) ||  handles.alf_SLR<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%-------------------------------------------------------------------------
function k_slr_Callback(hObject, eventdata, handles)
     handles.K_SLR=str2double(get(hObject,'String'));
    if isnan( handles.K_SLR) ||  handles.K_SLR<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%-------------------------------------------------------------------------
function n_plastic_slr_Callback(hObject, eventdata, handles)

   handles.n_plastic_SLR=str2double(get(hObject,'String'));
    n=handles.n_plastic_SLR;
    if isnan(n) || floor (n)~=n || n<=0
        errordlg(' Invalid input ! Input must be an integer.','Error','modal')
        set(hObject,'string','')
    end
    
    for i=1:n
        a{i}=strcat('p', num2str(i));
        b{1,i}='';
        b{2,i}='';
    end
    
    set (handles.table_slr, 'RowName', a)
    set (handles.table_slr, 'data', b')
    handles.table_SLR=str2double(b');
    handles.save_mat_index=0;
    guidata(hObject,handles)


function table_slr_CreateFcn(hObject, eventdata, handles)

function table_slr_CellEditCallback(hObject, eventdata, handles)
    if isnan(str2double(get(hObject,'data')))
        handles.table_SLR=(get(hObject,'data'));
    else
       handles.table_SLR=str2double(get(hObject,'data'));
    end
    guidata(hObject,handles)

%==========================================================================
% Shaft Material Properties
function n_plastic_shft_CreateFcn(hObject, eventdata, handles)
function k_shft_CreateFcn(hObject, eventdata, handles)
function alf_shft_CreateFcn(hObject, eventdata, handles)
function p_shft_CreateFcn(hObject, eventdata, handles)
function v_shft_CreateFcn(hObject, eventdata, handles)
function e_shft_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------------

function e_shft_Callback(hObject, eventdata, handles)
   handles.E_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.E_SHFT) || handles.E_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function v_shft_Callback(hObject, eventdata, handles)
    handles.v_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.v_SHFT) || handles.v_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')      
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function p_shft_Callback(hObject, eventdata, handles)
   handles.p_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.p_SHFT) || handles.p_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function alf_shft_Callback(hObject, eventdata, handles)
    handles.alf_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.alf_SHFT) || handles.alf_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function k_shft_Callback(hObject, eventdata, handles)
   handles.K_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.K_SHFT) || handles.K_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function n_plastic_shft_Callback(hObject, eventdata, handles)
    q=str2double(get(hObject,'String'));
    
    if isnan(q) || floor (q)~=q || q<=0
        errordlg(' Invalid input ! Input must be an integer.','Error','modal')
        set(hObject,'string','')
    else
       handles.n_plastic_SHFT=q; 
       n=q;
    end
    
    for i=1:n
        a{i}=strcat('p', num2str(i));
        b{1,i}='';
        b{2,i}='';
    end
    
    set (handles.table_shft, 'RowName', a)
    set (handles.table_shft, 'data', b')
    handles.table_SHFT=str2double(b');
    handles.save_mat_index=0;
    guidata(hObject,handles)
   
%--------------------------------------------------------------------------
function table_shft_CellEditCallback(hObject, eventdata, handles)
    
    if isnan(str2double(get(hObject,'data')))
        handles.table_SHFT=(get(hObject,'data'));
    else
        handles.table_SHFT=str2double(get(hObject,'data'));
    end
    
    handles.save_mat_index=0;
    guidata(hObject,handles)

%==========================================================================
function e_ins_CreateFcn(hObject, eventdata, handles)
function v_ins_CreateFcn(hObject, eventdata, handles)
function p_ins_CreateFcn(hObject, eventdata, handles)
function alf_ins_CreateFcn(hObject, eventdata, handles)
function k_ins_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------------

function e_ins_Callback(hObject, eventdata, handles)
   handles.E_INS=str2double(get(hObject,'String'));
    if isnan(handles.E_INS) || handles.E_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------

function v_ins_Callback(hObject, eventdata, handles)
   handles.v_INS=str2double(get(hObject,'String'));
    if isnan(handles.v_INS) || handles.v_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function p_ins_Callback(hObject, eventdata, handles)
    handles.p_INS=str2double(get(hObject,'String'));
    if isnan(handles.p_INS) || handles.p_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function alf_ins_Callback(hObject, eventdata, handles)
   handles.alf_INS=str2double(get(hObject,'String'));
    if isnan(handles.alf_INS) || handles.alf_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','') 
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------
function k_ins_Callback(hObject, eventdata, handles)
    handles.K_INS=str2double(get(hObject,'String'));
    if isnan(handles.K_INS) || handles.K_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')  
    end
    handles.save_mat_index=0;
    guidata(hObject,handles)
%--------------------------------------------------------------------------

function save_material_Callback(hObject, eventdata, handles)
    handles.table_SHFT
    b=sum(sum(isnan(handles.table_SHFT)));
    a=sum(sum(isnan(handles.table_SLR)));

    sl=0;
    sh=0;
    in=0;
    
    if isnan(handles.K_INS)||isnan(handles.alf_INS)||isnan(handles.p_INS) ||isnan(handles.v_INS)||isnan(handles.E_INS)
        in=1;
    end
    if isnan(handles.K_SLR)||isnan(handles.alf_SLR)||isnan(handles.p_SLR)||isnan(handles.v_SLR)||isnan(handles.E_SLR)|| a~=0
        sl=1;
    end
    if isnan(handles.K_SHFT)||isnan(handles.alf_SHFT)||isnan(handles.p_SHFT)||isnan(handles.v_SHFT)||isnan(handles.E_SHFT) || b~=0
        sh=1;
    end
    %---------------------------------
    str='';
    if sl==1
        str=strcat(str, ' "Slip Ring"');
    end
   
    if sh==1
        str=strcat(str, '  "Sahft" ');
    end
  
    if in==1
        str=strcat(str, '  "Insulation"');
    end
    %---------------------------------
    if  in==1 || sh==1 || sl==1
        if (in+sh+sl>1)
            str= strcat((strcat('Material properties of ', str)),' are incomplete!');
        else
             str= strcat((strcat('Material properties of ', str)),' is incomplete!');
        end
        errordlg(str,'Error','modal')  
    else
       msgbox('Data saved successfully !','Message','modal')
       handles.save_mat_index=1;   
    end   
 guidata(hObject,handles)   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Dimensions
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function dimen_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------------
function dimen_CellEditCallback(hObject, eventdata, handles) 
    handles.Dimen1=(get(hObject, 'data'));
    if isnan(str2double(handles.Dimen1(:,1)))
        handles.Dimen2=handles.Dimen1(:,1);
    else
       handles.Dimen2=str2double(handles.Dimen1(:,1));  
    end
    handles.save_dim_index=0;
    guidata(hObject, handles);
%--------------------------------------------------------------------------      
function save_dim_Callback(hObject, eventdata, handles)
    if sum(isnan(handles.Dimen2))>0
        errordlg('Incomplete data !','Error','modal')
    elseif  sum(handles.Dimen2<=0)>0
        errordlg('Invalid input! Dimensions Should be non-negative and non-zero.','Error','modal')    
    elseif floor(handles.Dimen2(end))~=handles.Dimen2(end) || floor(handles.Dimen2(end-1))~=handles.Dimen2(end-1)
        errordlg('Invalid input! Number of teeth and holes should be integer','Error','modal')  
    else
        msgbox('Data saved successfully !','Message','modal')
        handles.save_dim_index=1;
        
        n=handles.Dimen2(end);
            if n~=handles.n_old
                for i=1:n*4+1
                    T{i}=strcat('T', num2str(i-1));
                    D{i,1}='';
                    D{i,2}='';
                end
                set(handles. t_teeth, 'RowName', T);
                set(handles. t_teeth, 'data', D);
                handles.table_teeth2=str2double(D(:,1));
                handles.save_temp_index=0;
            end
    handles.n_old=handles.Dimen2(end);        
    end
    
guidata(hObject, handles);
%------------------
function help_dim_Callback(hObject, eventdata, handles)
close(figure(1));
j=figure(1);
set(j, 'Name', 'Dimensions')
set(j,'Position',[200 100 900 500],'Color', [.94 .94 .94])
ax1=axis()
 imshow('x.jpg')

    

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Temperature
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function theta_CreateFcn(hObject, eventdata, handles) 
function t_teeth_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------------
function t_teeth_CellEditCallback(hObject, eventdata, handles)
    handles.table_teeth1=(get(hObject,'data'));
    if isnan(str2double(handles.table_teeth1(:,1)))
        handles.table_teeth2=handles.table_teeth1(:,1);
    else
        handles.table_teeth2=str2double(handles.table_teeth1(:,1));
    end
    handles.save_temp_index=0;
    guidata(hObject,handles);
%--------------------------------------------------------------------------
function theta_CellEditCallback(hObject, eventdata, handles)
    handles.table_temp1=(get(hObject,'data'));
    
    if isnan(str2double(handles.table_temp1(:,1)))
        handles.table_temp2=handles.table_temp1(:,1);
    else
        handles.table_temp2=str2double(handles.table_temp1(:,1));
    end
    handles.save_temp_index=0;
   
    guidata(hObject,handles);
%--------------------------------------------------------------------------
function save_temp_Callback(hObject, eventdata, handles)
handles.table_teeth2;
handles.table_temp2;
    if  handles.save_dim_index==0
      errordlg('Please first complete dimensions and save it!','Error','modal')  
    elseif sum(isnan(handles.table_teeth2))~=0 || sum(isnan(handles.table_temp2))~=0
      errordlg('Data is incomplete!','Error','modal')  
    else
      handles.save_temp_index=1;  
      msgbox('Data saved successfully !','Message','modal')
    end
guidata(hObject,handles)
function help_temp_Callback(hObject, eventdata, handles)
close(figure(2));
j=figure(2);
set(j, 'Name', 'Temperature')

set(j,'Position',[200 100 700 500],'Color', [.94 .94 .94])
ax1=axis();

imshow('y.jpg');

    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Other Properties
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function fric_cof_CreateFcn(hObject, eventdata, handles)
function mesh_shft_CreateFcn(hObject, eventdata, handles)
function mass_conec_CreateFcn(hObject, eventdata, handles)
function bolt_pre_load_CreateFcn(hObject, eventdata, handles)
function rot_speed_CreateFcn(hObject, eventdata, handles)
function mesh_slr_CreateFcn(hObject, eventdata, handles)
function mesh_ins_CreateFcn(hObject, eventdata, handles)
function thermal_press_point_CreateFcn(hObject, eventdata, handles)
function i11_CreateFcn(hObject, eventdata, handles)
function i22_CreateFcn(hObject, eventdata, handles)
function i33_CreateFcn(hObject, eventdata, handles)

%-------------------------------------------------------------------------

function i11_Callback(hObject, eventdata, handles)
    handles.I11=str2double(get(hObject,'String'));
    if isnan(handles.I11) || handles.I11<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_other_index=0;
    guidata(hObject,handles);

%---------------------------------------------------------
function i22_Callback(hObject, eventdata, handles)
    handles.I22=str2double(get(hObject,'String'));
    if isnan(handles.I22) || handles.I22<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%---------------------------------------------------------
function i33_Callback(hObject, eventdata, handles)
    handles.I33=str2double(get(hObject,'String'));
    if isnan(handles.I33) || handles.I33<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%---------------------------------------------------------


function mesh_shft_Callback(hObject, eventdata, handles)
    handles.mesh_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.mesh_SHFT) || handles.mesh_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','')
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function mesh_slr_Callback(hObject, eventdata, handles)
    handles.mesh_SLR=str2double(get(hObject,'String'));
    if isnan(handles.mesh_SLR) || handles.mesh_SLR<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string',''); 
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function mesh_ins_Callback(hObject, eventdata, handles)
    handles.mesh_INS=str2double(get(hObject,'String'));
    if isnan(handles.mesh_INS) || handles.mesh_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal');
        set(hObject,'string',''); 
        
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function rot_speed_Callback(hObject, eventdata, handles)
   q=str2double(get(hObject,'String'));
   handles.rot_SPEED=q*2*3.14/60;  
    if isnan(q) || q<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal');
        set(hObject,'string','');
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function mass_conec_Callback(hObject, eventdata, handles)
    handles.mass_CONEC=str2double(get(hObject,'String'));
    if isnan(handles.mass_CONEC) || handles.mass_CONEC<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function bolt_pre_load_Callback(hObject, eventdata, handles)
   handles.bolt_LOAD=str2double(get(hObject,'String'));
    if isnan( handles.bolt_LOAD) ||  handles.bolt_LOAD<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string',''); 
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function fric_cof_Callback(hObject, eventdata, handles)
   handles.fric_COF=str2double(get(hObject,'String'));
    if isnan( handles.fric_COF) ||  handles.fric_COF<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function thermal_press_point_Callback(hObject, eventdata, handles)
    q=str2double(get(hObject,'String'));
    handles.therm_press_PNT=q; 
    n=q;
    if isnan(q) || floor (q)~=q || q<=0
        errordlg(' Invalid input ! Input must be an integer.','Error','modal')
            set(handles.thermal_press_point,'string','')
    end

    for i=1:n
        T{i}=strcat('p', num2str(i));
        D{i,1}='';
        D{i,2}='';
    end
    set(handles.thermal_pressure_table, 'RowName', T);
    set(handles.thermal_pressure_table, 'data', (D));
     handles.table_therm_press=str2double(D');
    handles.save_other_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function thermal_pressure_table_CellEditCallback(hObject, eventdata, handles)
    
    if isnan(str2double(get(hObject,'data')))
         handles.table_therm_press=(get(hObject,'data'));
    else
        handles.table_therm_press=str2double(get(hObject,'data'));
    end
    handles.save_other_index=0;
    guidata(hObject,handles);
    
%-------------------------------------------------------------------------- 
function save_other_prop_Callback(hObject, eventdata, handles)
    a=sum(sum(isnan(handles.table_therm_press)));
  
if  handles.save_dim_index==0
        errordlg('Please first complete dimensions and save it!','Error','modal')  

elseif isnan(handles.therm_press_PNT) ||...
        isnan(handles.fric_COF) ||isnan(handles.rot_SPEED) ||...
        isnan(handles.mesh_INS) ||isnan(handles.mesh_SLR) ||...
        isnan(handles.mesh_SHFT) || isnan(handles.bolt_LOAD) || isnan(handles.mass_CONEC)|| a>0

        errordlg('Data is incomplete!','Error','modal') 
else
        R1 = handles.Dimen2(18);
        R2 = handles.Dimen2(19);
        R3 =  handles.Dimen2(20);
        L3 =  handles.Dimen2(3);
        L5 = handles.Dimen2(5);
        L10 = handles.Dimen2(9);
        Rx=R3-(1.0*L5/L3)*(R3-R2);

        d1=floor((Rx-R1)*.95);
        d2=floor((L10)*.95);
        d3=min(d1,d2);

        if (d3<handles.mesh_INS)
            promptMessage = sprintf('The assigned Insulation mesh size may lead to a nonconformal mesh in the interfaces.Do you want to assigne an automatic value? ');
            titleBarCaption = 'Insulation mesh size';
            button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No','modal'); 

            if strcmp(button, 'Yes')
                set(handles.mesh_ins, 'string', num2str(d3))
                handles.mesh_INS=d3;
            end
        end

        handles.save_other_index=1;  
        msgbox('Data saved successfully !','Message','modal')

end
guidata(hObject,handles)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Fatigue S-N Data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function su_slr_CreateFcn(hObject, eventdata, handles)
function su_shft_CreateFcn(hObject, eventdata, handles)
function su_ins_CreateFcn(hObject, eventdata, handles)
function s_n_data_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------------- 
function su_slr_Callback(hObject, eventdata, handles)
    handles.su_SLR=str2double(get(hObject,'String'));
    if isnan(handles.su_SLR) || handles.su_SLR<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    handles.save_sn_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function su_shft_Callback(hObject, eventdata, handles)
    handles.su_SHFT=str2double(get(hObject,'String'));
    if isnan(handles.su_SHFT) || handles.su_SHFT<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    handles.save_sn_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function su_ins_Callback(hObject, eventdata, handles)
    handles.su_INS=str2double(get(hObject,'String'));
    if isnan(handles.su_INS) || handles.su_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    handles.save_sn_index=0;
    guidata(hObject,handles);

%-------------------------------------------------------------------------- 
function help_sn_Callback(hObject, eventdata, handles)
    close(figure(3));
    j=figure(3);
    set(j, 'Name', 'S-N Data')
    set(j,'Position',[200 100 900 500],'Color', [.94 .94 .94])
    ax1=axis();
    imshow('z.jpg');

%-------------------------------------------------------------------------- 
function s_n_data_CellEditCallback(hObject, eventdata, handles)
   
    if isnan(str2double(get(hObject,'data')))
         handles.S_N_Data=(get(hObject,'data'));
    else
         handles.S_N_Data=str2double(get(hObject,'data'));
    end
    
    handles.save_sn_index=0;
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function save_sn_Callback(hObject, eventdata, handles)
    if sum(sum(isnan(handles.S_N_Data)))>0 || isnan(handles.su_SLR) || isnan(handles.su_SHFT) || isnan(handles.su_INS)
        errordlg('Data is incomplete!','Error','modal') 

    else
        handles.save_sn_index=1;  
        msgbox('Data saved successfully !','Message','modal')
    end

guidata(hObject,handles)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Save and Load
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function file_Callback(hObject, eventdata, handles)

function save_Callback(hObject, eventdata, handles)
    uisave('handles','');
% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
   [filename1,filepath1]=uigetfile({'*.mat','.mat'},'Select Data File 1');
    cd(filepath1);
    handles0=load(filename1);
    handles1=handles0.handles;

    handles.E_SLR=handles1.E_SLR;              
    handles.v_SLR=handles1.v_SLR;             
    handles.p_SLR=handles1.p_SLR;
    handles.alf_SLR=handles1.alf_SLR;
    handles.K_SLR=handles1.K_SLR;
    handles.n_plastic_SLR=handles1.n_plastic_SLR;
    handles.table_SLR= handles1.table_SLR;
 
    handles.E_SHFT= handles1.E_SHFT;
    handles.v_SHFT=handles1.v_SHFT;
    handles.p_SHFT= handles1.p_SHFT;
    handles.alf_SHFT=handles1.alf_SHFT;
    handles.K_SHFT=handles1.K_SHFT;
    handles.n_plastic_SHFT= handles1.n_plastic_SHFT;
    handles.table_SHFT= handles1.table_SHFT;

    handles.E_INS= handles1.E_INS;
    handles.v_INS=handles1.v_INS;
    handles.p_INS= handles1.p_INS;
    handles.alf_INS=handles1.alf_INS;
    handles.K_INS=handles1.K_INS;
   
    handles.save_mat_index=handles1.save_mat_index;

    set(handles.e_slr,'string',handles.E_SLR);
    set(handles.v_slr,'string', handles.v_SLR);
    set(handles.p_slr,'string',handles.p_SLR);
    set(handles.k_slr,'string',handles.K_SLR);
    set(handles.alf_slr,'string',handles.alf_SLR);
    set(handles.n_plastic_slr,'string',handles.n_plastic_SLR);

    set(handles.e_shft,'string', handles.E_SHFT);
    set(handles.v_shft,'string', handles.v_SHFT);
    set(handles.p_shft,'string', handles.p_SLR);
    set(handles.k_shft,'string', handles.K_SHFT);
    set(handles.alf_shft,'string',handles.alf_SHFT);
    set(handles.n_plastic_shft,'string',handles.n_plastic_SHFT);

    set(handles.e_ins,'string',handles.E_INS);
    set(handles.v_ins,'string', handles.v_INS);
    set(handles.p_ins,'string', handles.p_INS);
    set(handles.k_ins,'string', handles.K_INS);
    set(handles.alf_ins,'string',handles.alf_INS);

    set(handles.table_slr, 'data', handles.table_SLR);
    set(handles.table_shft, 'data', handles.table_SHFT);
    

    %Dimension------------------------------------
    handles.Dimen1=handles1.Dimen1;
    handles.Dimen2=handles1.Dimen2;
    set(handles.dimen,'data', handles1.Dimen1);
    handles.save_dim_index=handles1.save_dim_index;
    handles.n_old=handles1.n_old;
    %temperature----------------------------------

    handles.save_temp_index=handles1.save_temp_index;
    handles.table_temp1=handles1.table_temp1;
    handles.table_teeth1= handles1.table_teeth1;
    handles.table_temp2=handles1.table_temp2;
    handles.table_teeth2= handles1.table_teeth2;
    
    set(handles.theta,'data', handles1.table_temp1);
    n=handles.Dimen2(end);
    for i=1:n*4+1
        T{i}=strcat('T', num2str(i-1));
    end
    set(handles. t_teeth, 'RowName', T);
    set(handles.t_teeth,'data', handles1.table_teeth1);
  
    % other data--------------------------------
    handles.save_other_index=handles1.save_other_index;
    handles.table_therm_press=handles1.table_therm_press;
    handles.therm_press_PNT=handles1.therm_press_PNT;
    handles.fric_COF=handles1.fric_COF;
    handles.bolt_LOAD=handles1.bolt_LOAD;
    handles.mass_CONEC=handles1.mass_CONEC;
    handles.rot_SPEED=handles1.rot_SPEED;
    handles.mesh_INS=handles1.mesh_INS;
    handles.mesh_SLR=handles1.mesh_SLR;
    handles.mesh_SHFT=handles1.mesh_SHFT;


    set(handles.fric_cof,'string',handles.fric_COF);
    set(handles.mass_conec,'string',handles.mass_CONEC);
    set(handles.bolt_pre_load,'string',handles.bolt_LOAD);
    set(handles.rot_speed,'string',handles.rot_SPEED*60/2/3.14);
    set(handles.mesh_slr,'string',handles.mesh_SLR);
    set(handles.mesh_ins,'string',handles.mesh_INS);
    set(handles.mesh_shft,'string',handles.mesh_SHFT);
    set(handles.thermal_press_point,'string',handles.therm_press_PNT);
    set(handles.thermal_pressure_table, 'data',  handles.table_therm_press);
    

    %SN Data--------------------------------------
    handles.S_N_Data=handles1.S_N_Data;
    set(handles.s_n_data,'data', handles1.S_N_Data);
    handles.save_sn_index=handles1.save_sn_index;
    handles.su_SLR=handles1.su_SLR;
    handles.su_SHFT=handles1.su_SHFT;
    handles.su_INS=handles1.su_INS;

    set(handles.su_slr,'string',handles.su_SLR);
    set(handles.su_shft,'string',handles.su_SHFT);
    set(handles.su_ins,'string', handles.su_INS);
    %Run------------------------------------------
    handles.low_cycle_Fatigue=handles1.low_cycle_Fatigue;
    handles.high_cycle_Fatigue= handles1.high_cycle_Fatigue;
    handles.mesh_Conv=handles1.mesh_Conv;
    handles.sep_Check=handles1.sep_Check;
    
    set(handles.low_cycle_fatigue, 'value', handles.low_cycle_Fatigue);
    set(handles.high_cycle_fatigue, 'value', handles.high_cycle_Fatigue);
    set(handles.mesh_conv, 'value', handles.mesh_Conv);
    set(handles.sep_check, 'value', handles.sep_Check);
    
    if handles.mesh_Conv==1
        set(handles.mesh_ref_fac,'Enable', 'on')
    else
        set(handles.mesh_ref_fac,'Enable', 'off')
    end
    guidata(hObject, handles);

    if handles.sep_Check==1
        set(handles.over_speed,'Enable', 'on')
    else
        set(handles.over_speed,'Enable', 'off')
    end

    if  handles.high_cycle_Fatigue==1
        set(handles.speed_tol,'Enable', 'on')
        set(handles.high_cycle_no,    'Enable', 'on')
    else
        set(handles.speed_tol,'Enable', 'off')
        set(handles.high_cycle_no,    'Enable', 'off')
    end
    
    if  handles.low_cycle_Fatigue==1
        set(handles.low_cycle_no,    'Enable', 'on')
    else
        set(handles.low_cycle_no,    'Enable', 'off')
    end
  
%   handles.model_Name=handles1.model_Name;
%   handles.save_Dir=handles1.save_Dir;
   handles.numCPU= handles1.numCPU;
   set(handles.num_cpu, 'string', handles.numCPU);
   handles.singular=handles1.singular;
     
    handles.LIN1=handles1.LIN1;
    handles.LIN2=handles1.LIN2;
    handles.Q1=handles1.Q1;
    handles.Q2=handles1.Q2;
    handles.Q3=handles1.Q3;
    set(handles.lin1,    'string', handles.LIN1)
    set(handles.lin2,    'string', handles.LIN2)
    set(handles.q1,    'string', handles.Q1)
    set(handles.q2,    'string', handles.Q2)
    set(handles.q3,    'string', handles.Q3) 

    handles.I11=handles1.I11;
    handles.I22=handles1.I22;
    handles.I33=handles1.I33;
    set(handles.i11,'string',handles.I11);
    set(handles.i22,'string',handles.I22);
    set(handles.i33,'string',handles.I33);
    
    guidata(hObject,handles)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Run
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function model_name_CreateFcn(hObject, eventdata, handles)
function save_dir_CreateFcn(hObject, eventdata, handles)
function mesh_ref_fac_CreateFcn(hObject, eventdata, handles)
function speed_tol_CreateFcn(hObject, eventdata, handles)
function over_speed_CreateFcn(hObject, eventdata, handles)
function high_cycle_no_CreateFcn(hObject, eventdata, handles)
function low_cycle_no_CreateFcn(hObject, eventdata, handles)
function num_cpu_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------------
function mesh_conv_Callback(hObject, eventdata, handles)
    handles.mesh_Conv=get(hObject, 'value');
    
    if handles.mesh_Conv==1
        set(handles.mesh_ref_fac,'Enable', 'on')
    else
        set(handles.mesh_ref_fac,'Enable', 'off')
    end
    guidata(hObject, handles);
%-------------------------------------------------------------------------
function sep_check_Callback(hObject, eventdata, handles)
    handles.sep_Check=get(hObject, 'value');
    
    if handles.sep_Check==1
        set(handles.over_speed,'Enable', 'on')
    else
        set(handles.over_speed,'Enable', 'off')
    end
    guidata(hObject, handles);
%-------------------------------------------------------------------------    
function high_cycle_fatigue_Callback(hObject, eventdata, handles)
   handles.high_cycle_Fatigue=get(hObject, 'value');
    
    if  handles.high_cycle_Fatigue==1
        set(handles.speed_tol,'Enable', 'on')
        set(handles.high_cycle_no,    'Enable', 'on')
       
    else
        set(handles.speed_tol,'Enable', 'off')
        set(handles.high_cycle_no,    'Enable', 'off')
    end
    guidata(hObject, handles);
%-------------------------------------------------------------------------    
function low_cycle_fatigue_Callback(hObject, eventdata, handles)
   handles.low_cycle_Fatigue=get(hObject, 'value');
    
    if  handles.low_cycle_Fatigue==1
        
        set(handles.low_cycle_no,    'Enable', 'on')
       
    else
        
        set(handles.low_cycle_no,    'Enable', 'off')
    end
    guidata(hObject, handles);
    

%-------------------------------------------------------------------------

function brows_save_dir_Callback(hObject, eventdata, handles)
    [path] = uigetdir();
    path(path=='\')='/';
    handles.save_Dir=path;
    if (path)
        set(handles.save_dir, 'string', path);
        handles.save_Dir=path;

    else
        errordlg('Please Select a Directory!','Error','modal') 
    end
    guidata(hObject, handles)

%--------------------------------------------------------------------------
function speed_tol_Callback(hObject, eventdata, handles)
    q=str2double(get(hObject, 'string'));
    handles.speed_Tol=q*2*3.14/60;
    guidata(hObject, handles);
%--------------------------------------------------------------------------

function mesh_ref_fac_Callback(hObject, eventdata, handles)
    handles.mesh_ref_Fac=str2double(get(hObject, 'string'));
    guidata(hObject, handles);
%--------------------------------------------------------------------------
function over_speed_Callback(hObject, eventdata, handles)
q=str2double(get(hObject, 'string'));
handles.over_Speed=q*2*3.14/60;
guidata(hObject, handles);
%--------------------------------------------------------------------------     
function high_cycle_no_Callback(hObject, eventdata, handles)
    handles.high_cycle_No=str2double(get(hObject,'string'));
    
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function low_cycle_no_Callback(hObject, eventdata, handles)
    handles.low_cycle_No=str2double(get(hObject,'string'));
  
    guidata(hObject,handles);
    
%--------------------------------------------------------------------------
function save_dir_Callback(hObject, eventdata, handles)


%--------------------------------------------------------------------------

function model_name_Callback(hObject, eventdata, handles)
    handles.model_Name=get(hObject, 'string');
    guidata(hObject, handles);
    
function num_cpu_Callback(hObject, eventdata, handles)
    handles.numCPU=str2double(get(hObject,'String'));
    n=handles.numCPU;
    if isnan(n) || floor (n)~=n || n<=0
        errordlg(' Invalid input ! Input must be an integer.','Error','modal')
        set(hObject,'string','')
    end
guidata(hObject,handles);    

function q1_CreateFcn(hObject, eventdata, handles)
function q3_CreateFcn(hObject, eventdata, handles)
function q2_CreateFcn(hObject, eventdata, handles)
function lin1_CreateFcn(hObject, eventdata, handles)
function lin2_CreateFcn(hObject, eventdata, handles)
%----------------------------------------------------------------------------
function q1_Callback(hObject, eventdata, handles)
    handles.Q1=str2double(get(hObject,'String'));
    if isnan(handles.su_INS) || handles.su_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    guidata(hObject,handles);
%----------------------------------------------------------------------------
function q3_Callback(hObject, eventdata, handles)
    handles.Q2=str2double(get(hObject,'String'));
    if isnan(handles.su_INS) || handles.su_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    guidata(hObject,handles);
%----------------------------------------------------------------------------
function q2_Callback(hObject, eventdata, handles)
    handles.Q1=str2double(get(hObject,'String'));
    if isnan(handles.su_INS) || handles.su_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    guidata(hObject,handles);
%----------------------------------------------------------------------------
function lin1_Callback(hObject, eventdata, handles)
    handles.LIN1=str2double(get(hObject,'String'));
    if isnan(handles.su_INS) || handles.su_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    guidata(hObject,handles);
%----------------------------------------------------------------------------
function lin2_Callback(hObject, eventdata, handles)

    handles.LIN2=str2double(get(hObject,'String'));
    if isnan(handles.su_INS) || handles.su_INS<=0
        errordlg('Invalid input! Data Should be non-negative and non-zero double.','Error','modal')
        set(hObject,'string','');
    end
    guidata(hObject,handles);

%--------------------------------------------------------------------------

function run_Callback(hObject, eventdata, handles)
promptMessage = sprintf('Select "Yes" for Analysis in Abaqus GUI, "No" for Analysis without Abaqus GUI. ');
titleBarCaption = 'Analysis Mode';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Cancel','modal');
if strcmp(button, 'Yes') || strcmp(button, 'No')
    
    incomplete={'Folowing Data is Incomplete:' };
    id=2;
     sn=0; mat=0; dim=0; temp=0; other=0; hcf=0; lcf=0; sepc=0; dirc=0; mna=0; sing=0;
    
 handles.save_mat_index;
    
    if  handles.save_mat_index==0
       mat=1;
       incomplete{id}='*Material Properties'
       id=id+1;
    end
    if handles.save_dim_index==0
       dim=1;
       incomplete{id}='*Dimensions'
       id=id+1;
    end
    if handles.save_temp_index==0
       temp=1
       incomplete{id}='*Temperatures'
       id=id+1;
    end
    
    if handles.save_other_index==0
       other=1;
       incomplete{id}='*Other Properties'
       id=id+1;
    end
    
    if handles.save_sn_index==0   && ( handles.high_cycle_Fatigue==1||handles.low_cycle_Fatigue==1)
       sn=1;
       incomplete{id}='*S-N-Data'
       id=id+1;
    end
    
    if  handles.high_cycle_Fatigue==1 &&(isnan(handles.high_cycle_No)|| isnan(handles.speed_Tol))
       hcf=1;
       incomplete{id}='*High Cycle No. or Speed Tol'
       id=id+1;
    end
    
    if  handles.low_cycle_Fatigue==1 &&(isnan(handles.low_cycle_No))
       lcf=1;
       incomplete{id}='*Low Cycle No.'
       id=id+1;
    end
    
    if  handles.sep_Check==1 &&(isnan(handles.over_Speed))
       sepc=1;
       incomplete{id}='*Over Speed'
       id=id+1;
    end

    if  handles.mesh_Conv==1 &&(isnan(handles.mesh_ref_Fac))
       sepc=1;
       incomplete{id}='*Mesh Refinment Factor'
       id=id+1;
    end

    if  strcmp(handles.save_Dir, '')
       dirc=1;
       incomplete{id}='*Save Dir.'
       id=id+1;
    end

    if  strcmp(handles.model_Name, '')
       mna=1;
       incomplete{id}='*Model Name'
       id=id+1;
    end    

%     if  isnan(handles.Q1)|| isnan(handles.Q2)||isnan(handles.Q3)||isnan(handles.lin1)||isnan(handles.lin2)
%        sing=1;
%        incomplete{id}='*Quad Distances or Lin. Distances'
%     end   

    if  (sn+other+dim+mat+temp+lcf+hcf+sepc+dirc+mna+sing)>0  
        errordlg(incomplete,'Error','modal')
    else
        
        handles.DIRODB=strcat(handles.save_Dir,'/' ,handles.model_Name,'Results');
        handles.DIRCAE=strcat(handles.save_Dir,'/' ,handles.model_Name,'MODELS');
        mkdir( handles.DIRODB)
        mkdir( handles.DIRCAE)
        
        fid = fopen('Var.py', 'w');
    %     fprintf(fid,'run_type = %d\n',handles.run_Type);
        fprintf(fid,'L1 = %0.12f\n',  handles.Dimen2(1));
        fprintf(fid,'L2 = %0.12f\n' , handles.Dimen2(2));
        fprintf(fid,'L3 = %0.12f\n',  handles.Dimen2(3));
        fprintf(fid,'L4 = %0.12f\n',  handles.Dimen2(4));
        fprintf(fid,'L5 = %0.12f\n',  handles.Dimen2(5));
        fprintf(fid,'L6 = %0.12f\n',  handles.Dimen2(6));
        fprintf(fid,'L8 = %0.12f\n',  handles.Dimen2(7));
        fprintf(fid,'L9 = %0.12f\n',  handles.Dimen2(8));
        fprintf(fid,'L10 = %0.12f\n', handles.Dimen2(9));
        fprintf(fid,'L11 = %0.12f\n', handles.Dimen2(10));
        fprintf(fid,'L12 = %0.12f\n', handles.Dimen2(11));
        fprintf(fid,'L16 = %0.12f\n', handles.Dimen2(12));
        fprintf(fid,'L17 = %0.12f\n', handles.Dimen2(13));
        fprintf(fid,'L18 = %0.12f\n', handles.Dimen2(14));
        fprintf(fid,'L19 = %0.12f\n', handles.Dimen2(15));
        fprintf(fid,'L20 = %0.12f\n', handles.Dimen2(16));
        fprintf(fid,'L21 = %0.12f\n', handles.Dimen2(17));
        fprintf(fid,'R1 = %0.12f\n',  handles.Dimen2(18));
        fprintf(fid,'R2 = %0.12f\n',  handles.Dimen2(19));
        fprintf(fid,'R3 = %0.12f\n',  handles.Dimen2(20));
        fprintf(fid,'R4 = %0.12f\n',  handles.Dimen2(21));
        fprintf(fid,'R5 = %0.12f\n',  handles.Dimen2(22));
        fprintf(fid,'R6 = %0.12f\n',  handles.Dimen2(23));
        fprintf(fid,'R7 = %0.12f\n',  handles.Dimen2(24));
        fprintf(fid,'R8 = %0.12f\n',  handles.Dimen2(25));
        fprintf(fid,'R9 = %0.12f\n',  handles.Dimen2(26));
        fprintf(fid,'R10 = %0.12f\n', handles.Dimen2(27));
        fprintf(fid,'D1 = %0.12f\n',  handles.Dimen2(28));
        fprintf(fid,'D2 = %0.12f\n',  handles.Dimen2(29));
        fprintf(fid,'D3 = %0.12f\n',  handles.Dimen2(30));
        fprintf(fid,'nH = %d\n',      handles.Dimen2(31));
        fprintf(fid,'n = %d\n',       handles.Dimen2(32));
        fprintf(fid,'#**********************************************\n');
        fprintf(fid,'Eslp = %0.12f\n',    handles.E_SLR);
        fprintf(fid,'vslp = %0.12f\n',    handles.v_SLR);
        fprintf(fid,'Exp_slp= %0.12f\n',  handles.alf_SLR);
        fprintf(fid,'Den_slp= %d\n',      handles.p_SLR);
        fprintf(fid,'k_slp = %d\n',       handles.K_SLR);
        fprintf(fid,'Eshaft = %0.12f\n',  handles.E_SHFT);
        fprintf(fid,'vshaft = %0.12f\n',  handles.v_SHFT);
        fprintf(fid,'Den_sh = %0.12f\n',  handles.p_SHFT);
        fprintf(fid,'Exp_sh= %d\n',       handles.alf_SHFT);
        fprintf(fid,'k_sh = %d\n',        handles.K_SHFT);  
        fprintf(fid,'Eins = %0.12f\n',    handles.E_INS);
        fprintf(fid,'vins = %0.12f\n',    handles.v_INS);
        fprintf(fid,'Den_ins= %d\n',      handles.p_INS);
        fprintf(fid,'Exp_ins = %d\n',     handles.alf_INS);   
        fprintf(fid,'k_ins = %d\n',       handles.K_INS);   
        
        fprintf(fid,'Syslp=[');
        for i=1: handles.n_plastic_SLR
            fprintf(fid,'[');
            fprintf(fid,'%f, %f ',handles.table_SLR(i,2), handles.table_SLR(i, 1));
            fprintf(fid,'],');
        end
        fprintf(fid,']\n');

        fprintf(fid,'Sysh=[');
        for i=1: handles.n_plastic_SHFT
            fprintf(fid,'[');
            fprintf(fid,'%f, %f ', handles.table_SHFT(i,2),  handles.table_SHFT(i, 1));
            fprintf(fid,'],');
        end
        fprintf(fid,']\n');
        fprintf(fid,'#**********************************************\n');
        fprintf(fid,'teta1 = %0.12f\n', handles.table_temp2(1));
        fprintf(fid,'teta2 = %0.12f\n', handles.table_temp2(2));
        fprintf(fid,'teta3 = %0.12f\n', handles.table_temp2(3));
        fprintf(fid,'teta4 = %0.12f\n', handles.table_temp2(4));
        fprintf(fid,'teta5 = %0.12f\n', handles.table_temp2(5));
        fprintf(fid,'teta6 = %0.12f\n', handles.table_temp2(6));
        fprintf(fid,'teta7 = %0.12f\n', handles.table_temp2(7));
        fprintf(fid,'teta8 = %0.12f\n', handles.table_temp2(8));
        fprintf(fid,'Tshaft = %0.12f\n',handles.table_temp2(9));   

        fprintf(fid,'T=[');
        for i=1: handles.Dimen2(end)*4+1
            fprintf(fid,'%f', handles.table_teeth2(i));
             fprintf(fid,', '); 
        end
        fprintf(fid,']\n');   
        fprintf(fid,'#**********************************************\n');

        fprintf(fid,'shaftMeshSize = %0.12f\n', handles.mesh_SHFT);
        fprintf(fid,'slipMeshSize = %0.12f\n', handles.mesh_SLR);
        fprintf(fid,'shrinkMeshSize = %0.12f\n',  handles.mesh_INS);
        fprintf(fid,'FC = %0.12f\n',handles.fric_COF);
        fprintf(fid,'Mass = %0.12f\n',handles.mass_CONEC);
        fprintf(fid,'w = %0.12f\n', handles.rot_SPEED);
        fprintf(fid,'boltLoad = %d\n',       handles.bolt_LOAD);
        fprintf(fid,'i11 = %d\n',       handles.I11);
        fprintf(fid,'i22 = %d\n',       handles.I22);
        fprintf(fid,'i33 = %d\n',       handles.I33);
        
        
        fprintf(fid,'prsscon=[');
        for i=1: handles.therm_press_PNT
            fprintf(fid,'[');
            fprintf(fid,'%f, %f ', handles.table_therm_press(i,2), handles.table_therm_press(i, 1));
            fprintf(fid,'],');
        end
        fprintf(fid,']\n');
        fprintf(fid,'#**********************************************\n');

        if handles.low_cycle_Fatigue==1|| handles.high_cycle_Fatigue==1
            fprintf(fid,'S_N_SLR= [');
            fprintf(fid,'(%0.1f,%0.1f, %0.1f, %0.1f),',handles.S_N_Data(1,1),handles.S_N_Data(1,3),handles.S_N_Data(1,5),handles.S_N_Data(1,7));
            fprintf(fid,'(%0.1f,%0.1f, %0.1f, %0.1f),',handles.S_N_Data(1,2),handles.S_N_Data(1,4),handles.S_N_Data(1,6),handles.S_N_Data(1,8));
            fprintf(fid,']\n');

            fprintf(fid,'S_N_SHFT= [');
            fprintf(fid,'(%0.1f,%0.1f, %0.1f, %0.1f),',handles.S_N_Data(2,1),handles.S_N_Data(2,3),handles.S_N_Data(2,5),handles.S_N_Data(2,7));
            fprintf(fid,'(%0.1f,%0.1f, %0.1f, %0.1f),',handles.S_N_Data(2,2),handles.S_N_Data(2,4),handles.S_N_Data(2,6),handles.S_N_Data(2,8));
            fprintf(fid,']\n');

            fprintf(fid,'S_N_INS= [');
            fprintf(fid,'(%0.1f,%0.1f, %0.1f, %0.1f),',handles.S_N_Data(3,1),handles.S_N_Data(3,3),handles.S_N_Data(3,5),handles.S_N_Data(3,7));
            fprintf(fid,'(%0.1f,%0.1f, %0.1f, %0.1f),',handles.S_N_Data(3,2),handles.S_N_Data(3,4),handles.S_N_Data(3,6),handles.S_N_Data(3,8));
            fprintf(fid,']\n');
             
            fprintf(fid,'S_U_SLR = %0.12f\n',  handles.su_SLR);
            fprintf(fid,'S_U_SHFT = %0.12f\n', handles.su_SHFT);
            fprintf(fid,'S_U_INS = %0.12f\n',  handles.su_INS);
        end



        if  handles.low_cycle_Fatigue ==1

            fprintf(fid,'N_LowCycle = %0.12f\n',handles.low_cycle_No);
        end

       if  handles.high_cycle_Fatigue ==1

           fprintf(fid,'N_HighCycle = %0.12f\n',handles.high_cycle_No);
           fprintf(fid,'w_tol = %0.12f\n',handles.speed_Tol);
       end
        if handles.mesh_Conv ==1

           fprintf(fid,'mesh_ref_fac = %0.12f\n',handles.mesh_ref_Fac);
        end
        if  handles.sep_Check ==1

           fprintf(fid,'w_overSpeed = %0.12f\n',handles.over_Speed);


        end

        fprintf(fid,'jobName= "%s" \n',handles.model_Name);
        fprintf(fid,'directoryname= "%s"\n', handles.DIRCAE);
        fprintf(fid,'directory= "%s"\n', handles.DIRODB);
        fprintf(fid,'numCPU= %d\n',handles.numCPU);
        fprintf(fid,'singpntCriteria= "%s"\n',handles.singular);
        
        fprintf(fid,'seperationCheck= %d\n',handles.sep_Check);
        fprintf(fid,'meshConvergence= %d\n',handles.mesh_Conv);
        fprintf(fid,'lowCycleAnalysis= %d\n',handles.low_cycle_Fatigue);
        fprintf(fid,'highCycleAnalysis= %d\n',handles.high_cycle_Fatigue);
         
        fprintf(fid,'quadDistance= [');
        fprintf(fid,'%0.1f, %0.1f, %0.1f',handles.Q1,handles.Q2,handles.Q3);
        fprintf(fid,']\n');
        
        fprintf(fid,'linDistance= [');
        fprintf(fid,'%0.1f, %0.1f',handles.LIN1,handles.LIN2);
        fprintf(fid,']\n');
        
        
        if  strcmp(button, 'No')
            system(['abaqus cae ','noGUI','=Functions1.py']);
        else 
            system(['abaqus cae ','script','=Functions1.py']);
        end
        
        pwd0=pwd;
        cd(handles.DIRODB);
        delete('*.ipm','*.prt','*.sim','*.ipm','*.sta','*.log','*.inp','*.com','*.dat','*.msg','*.odb_f')
        delete('*.rpy','*.1','*.2','*.pac','*.pac','*.res','*.sel','*.stt','*.rec','*.rpy','*.log','*.loc','*.abq','*.mdl')
        cd( pwd0);
        delete('*.1', '*.2','*.log', '*.rpy', '*.asv')
    end
    
    Infolder=dir(handles.DIRODB);
    MyListOfFiles = {Infolder(~[Infolder.isdir]).name};
    set(handles.listbox5, 'string', MyListOfFiles);

    Infolder      = dir(handles.DIRCAE);
    MyListOfFiles = {Infolder(~[Infolder.isdir]).name};
    set(handles.listbox6, 'string', MyListOfFiles);

    set(handles.uipanel_run,'visible', 'off') 
    set(handles.uipanel_results,'visible', 'on') 
    
guidata(hObject,handles);
end




%results panel------------------------------------------------------------
function result_tog_Callback(hObject, eventdata, handles)
 
set(handles.uipanel_run,'visible', 'off') 
set(handles.uipanel_results,'visible', 'on') 
set(handles.result_tog,'value', 0) 



guidata(hObject,handles);
%-------------------------------------------------------------------------- 
function run_tog_Callback(hObject, eventdata, handles)
    set(handles.uipanel_run,'visible', 'on') 
    set(handles.uipanel_results,'visible', 'off')
    set(handles.run_tog,'value', 0) 
    guidata(hObject,handles);
%--------------------------------------------------------------------------  
function listbox6_CreateFcn(hObject, eventdata, handles)
function listbox5_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------------- 
function listbox5_Callback(hObject, eventdata, handles)
    winopen(handles.DIRODB);
    guidata(hObject,handles);

%-------------------------------------------------------------------------- 
function listbox6_Callback(hObject, eventdata, handles)
    winopen( handles.DIRCAE);
    guidata(hObject,handles);
%-------------------------------------------------------------------------- 
