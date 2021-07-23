VNI <- read_csv("VNI.csv", col_types = cols(`Change %` = col_skip(),
                                            High = col_skip(), Low = col_skip(),
                                            Open = col_skip(), Vol. = col_skip()))
Date = VNI$Date
dateNew <- vector()
for(i in 1:length(Date)){dateNew[i] = as.Date.character(paste(substr(Date[i],1,3),"/",substr(Date[i],5,6),"/",substr(Date[i],9,12),sep = ""),format = "%b/%d/%Y")}
vni = zoo::zoo(x = VNI$Price,order.by = as.Date.numeric(dateNew,origin = "1970-01-01"))
