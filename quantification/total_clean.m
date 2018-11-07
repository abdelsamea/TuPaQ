function clean=total_clean(bwnuclei,fill)

if fill == 1
    warning off;
    bwnuclei=imfill(bwnuclei,'holes');
    warning on;
end




clean=bwnuclei;
