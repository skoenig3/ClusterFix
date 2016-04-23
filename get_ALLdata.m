function [time_arr,event_arr,eog_arr,epp_arr, header,trialcount]  = get_ALLdata(name)

% function to read the cortex datafile 
% and put the timecodes, eventcodes, trialheaders
% in an array of arrays, and the nr of trials in a scalar.
% - - -
% USAGE: [times,events,eog_arr,epp_arr, header,trialcount] = get_data(path\file_name)
% - - - 
% adapted from readcort.m
%
% GDLH modified this 5/16/00 -- no need to preallocate the sizes of the arrays
% Matlab 5.2 will automatically pad with zeros when necessary.  How convenient!
% (But preallocating does seem to speed things up a bit.)
%
% GDLH modified this 4/2/01 -- Now returning the eog_arr as well.
%
% GDLH modified this 4/28/01 -- if we hit an error in the datafile, just
% return what we can.

MAXNUMTRIALS = 1000;
trialcount = 0;
fid = fopen(name, 'rb');
if (fid == -1)
   error(['Cannot open file: ',name]);
end
time_arr  = zeros(1, MAXNUMTRIALS);	% array of time_code arrays
event_arr = zeros(1, MAXNUMTRIALS);	% array of event_code arrays
eog_arr = zeros(1, MAXNUMTRIALS);	% array of eog_code arrays
eog_lengths = [];
epp_arr = zeros(1, MAXNUMTRIALS);	% array of epp_code arrays
epp_lengths = [];
header    = zeros(13, MAXNUMTRIALS);	% array of trial headers
hd=zeros(1,13);					% a single header (read every trial)
while (~feof (fid))
   length = fread(fid, 1, 'ushort');
   if (isempty(length)~=1)
      hd(1,1:8)= (fread(fid, 8, 'ushort'))'; 
      hd(1,9:10) = (fread(fid, 2, 'uchar'))'; 
      hd(1,11:13)= (fread(fid, 3, 'ushort'))';
      hd(1,5)=hd(1,5)/4;	% reduce the size of the time_code array 
      hd(1,6)=hd(1,6)/2;	% idem for event_code array
      hd(1,7)=hd(1,7)/2;  	% and eogs
      hd(1,8)=hd(1,8)/2;	% and epp (not used)
      % read the time codes, event_codes, eogs and epps (if any) for this trial
      try
         time_arr(1:(hd(1,5)),trialcount+1)   = (fread (fid,(hd(1,5)) , 'ulong'));   
         event_arr(1:(hd(1,6)),trialcount+1) = (fread (fid,(hd(1,6)), 'ushort'));
         
         % epp array is stored in the data file before the eog array
         epp_arr(1:(hd(1,8)),trialcount+1)        = fread (fid,(hd(1,8)), 'short');
         
         % must modify all epp entries
			for i = 1:1:(hd(1,8)),
        
        	% Must extract the data from the raw epp values.
         % The epp value is made up of 12-bits of data, and 4-bits (the
   	      % low-order 4 bits) of the channel number.  
	         % To extract the data, must right shift the raw data by 4 bits (to
   	      % get rid of the channel number and put the data value in the proper 
      	   % location).  After this conversion, you must still add or subtract
         	% 2048, since otherwise the value is not right.  (I think that this
	         % is because of the way that matlab handles negative values during
   	      % the bitwise operations.)
      	   % These calculations cause the results of ctx2txt.m to be the same as
         	% for cortview.exe for the EOG and EPP values.
         
         	eppdata = bitshift(epp_arr(i, trialcount+1), -4);
                  
	         if (eppdata < 0)
   	         eppdata = eppdata  + 2047;
      	   else
         	   eppdata = eppdata - 2048;
	         end;
         
            epp_arr(i, trialcount+1) = eppdata;
            
                      
                    
         end;
         epp_lengths(trialcount+1) = hd(1,8);

         
         eog_arr(1:(hd(1,7)),trialcount+1)     = (fread (fid,(hd(1,7)), 'short'));
         eog_lengths(trialcount+1) = hd(1,7);
                  
         % put this trial's header in the header_array
         header(1:13,trialcount+1)=hd(1,1:13)';
         trialcount = trialcount+1;
      catch
         break;
      end;
   end; 
end;

fclose(fid);
time_arr(:,[trialcount+1:end]) = [];
event_arr(:,[trialcount+1:end]) = []; 
eog_arr(:,[trialcount+1:end]) = []; 
epp_arr(:,[trialcount+1:end]) = []; 
header(:,[trialcount+1:end]) = []; 

if (size(eog_arr,1) == 1)
   eog_arr = [];
else
   % 'nan' padding the eog_ary
   for i = 1:trialcount
      %i
      %eog_lengths(i)
      %size(eog_arr)
      if (eog_lengths(i) < size(eog_arr,1))
         eog_arr([eog_lengths(i)+1:end],i) = nan;   
      end
   end
end

if (size(epp_arr,1) == 1)
   epp_arr = [];
else
   % 'nan' padding the epp_arr
   for i = 1:trialcount
      if (epp_lengths(i) < size(epp_arr,1))
        epp_arr([epp_lengths(i)+1:end],i) = nan;   
      end
   end
end

