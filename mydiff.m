function rslt = mydiff(data)
rslt = diff(data);
rslt = [rslt(1,:); rslt];
end
