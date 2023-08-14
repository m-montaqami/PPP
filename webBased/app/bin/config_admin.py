
Info={
    'company':'KNTU',
    'author':'m-montaqami',
}


software={
    'version':'1.0',
    'remote-mode':'offline',
    'OS':'linux'
}

Files={
    'EGM':"app/bin/geoids/egm2008-5.pgm",
}

Ex_Links={
    'RequestTimeOut':2,
    'itrg':'http://chapman.upc.es/irtg/last_results/MOST-RECENT-RT-SNAPSHOT-IONEX-FILE.Z',
    'upc':'http://chapman.upc.es/tomion/real-time/quick/last_results.uadg/MOST-RECENT-RT-DAILY-IONEX-FILE.Z',
    'ECMWF':{
    'product':'reanalysis-era5-pressure-levels',
    'request':{
        'product_type': 'reanalysis',
        'variable': ["geopotential","temperature"],
        'pressure_level': [
        '1', '2', '3','5', '7', '10','20', '30', '50','70', '100', '125','150', '175', '200',
        '225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750',
        '775', '800', '825','850', '875', '900','925', '950', '975','1000'],
        'year': '2012',
        'month': '12',
        'day': '12',
        'time':[],
        'area': [53, 33, 49,37], 
        'format': 'grib'
        },
    'downloadPath':'archive/trop/ECMWF/ecmwf.grib'
},
    'GFS':{
        'url':['https://nomads.ncep.noaa.gov/cgi-bin/filter_fnl.pl?file=gdas.t','z.pgrb2.1p00.f009&lev_surface=on&var_PRES=on&dir=%2Fgdas.',
                '%2F','%2Fatmos'],
        'item':['gdas.t','z.pgrb2.1p00.f009'],
        'downloadPath':'archive/trop/GFS/gfs.grib'
        
    }

}
