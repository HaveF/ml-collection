function resultlabel = resign_region_label(finres)
tempfinres = finres;
original_max_label = max(max(tempfinres));
resultlabel = zeros(size(tempfinres));
baselabel =0;
for i=1:original_max_label
    mregion = tempfinres==i;
    label = bwlabel(mregion);
    maxlabel = max(max(label));
    label = label + baselabel;
    label(find(tempfinres~=i))=0;
    baselabel = baselabel + maxlabel;
    resultlabel = resultlabel + label;    
end
end