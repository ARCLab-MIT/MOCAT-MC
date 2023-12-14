function satellite_type = filter_objclass_fragments_int(objclass_org)
%Set object class to one of these satellite types ['Payload Fragmentation Debris'=3,'Rocket Fragmentation Debris'=7,'Other Debris'=10 or 'Unknown'=11]
%after explosion or collision.
%objclass: 
%         1 = 'Payload'
%         2 = 'Payload Mission Related Object' 
%         3 = 'Payload Fragmentation Debris'
%         4 = 'Payload Debris' 
%         5 = 'Rocket Body' 
%         6 = 'Rocket Mission Related Object' 
%         7 = 'Rocket Fragmentation Debris' 
%         8 = 'Rocket Debris'
%         9 = 'Debris' (new) 
%         10 = 'Other Debris' 
%         11 = 'Unknown' or []
    
if objclass_org<4.5
    satellite_type = 3; %'Payload Fragmentation Debris'
elseif objclass_org<8.5
    satellite_type = 7; %'Rocket Fragmentation Debris' 
elseif objclass_org==9
    satellite_type = 9; %'Debris'
elseif objclass_org==10
    satellite_type = 10; %'Other Debris'
elseif objclass_org==11
    satellite_type = 11; %'Unknown' or []
else 
    objclass_org = objclass_org
    error('objclass did not match any of the pre-determined options. Please review this function to assign object class to one of these satellite types: [''Payload Fragmentation Debris'',''Rocket Fragmentation Debris'',''Debris'', ''Other Debris'' or ''Unknown'']')
end

