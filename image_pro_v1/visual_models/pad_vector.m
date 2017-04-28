function nv=pad_vector(v,PAD)

% Ensure continous vector
v(find(v==1,1,'first'):find(v==1,1,'last'))=1;

vd=find(diff(v)~=0);
if length(vd)>1
    if vd(1)-PAD<1, y1=1; else y1=vd(1)-PAD; end
    if vd(end)+PAD>length(v), y2=length(v); else y2=vd(end)+PAD; end  
elseif length(vd)==1
   if v(1)==1, 
       y1=1; 
   elseif vd(1)-PAD<1, 
       y1=1; 
   else
       y1=vd(1)-PAD;
   end
   
   if v(end)==1
       y2=length(v); 
   elseif vd(1)+PAD>length(v), 
       y2=length(v); 
   else
       y2=vd(1)+PAD;
   end
else
   y1=1; y2=length(v); 
end
nv=zeros(1,length(v));    
nv(y1:y2)=1;
nv=logical(nv);

end