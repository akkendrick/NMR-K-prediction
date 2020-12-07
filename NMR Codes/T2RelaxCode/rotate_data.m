function [data_rot,noise,n]=rotate_data(data)
if size(data,1)<=size(data,2)
    data=data';
end

 y1=data;
 
 if length(y1) < 30
     endIndex = length(y1);
 else
     endIndex = 30;
 end
    ne=length(y1);
    phi=atan2(sum(imag(y1(1:endIndex))),sum(real(y1(1:endIndex))));
    y1=(ones(ne(1),1)*exp(-1i*phi')).*y1;
    n=std(imag(y1),1);
    data_rot=real(y1);
    noise=imag(y1);

end