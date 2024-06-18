// Import 3 assets (from the anizam1 project) as inputs:
// an 'ecosystem' table that is a rectangular extent of the area of interest
// a 'gmas_vec' table of the GMA polygons
// a 'parks_vec' table of the national park polygons

// Code just below (until two blank comment lines) illustrates creation of these assets
// Coordinates for the bounds of a rectangle.
var xMin = 25.42;
var yMin = -13.81;
var xMax = 26.58;
var yMax = -15.54;

var gke_rect = ee.Geometry.Rectangle(
  [
    [xMin, yMin],
    [xMax, yMax]   // max x and y
  ]
);
Map.addLayer(gke_rect, {}, 'gke_rect');

var gke_extent = ee.FeatureCollection(gke_rect);

Export.table.toAsset({
  collection: gke_extent,
  description:'create_gke_extent',
  assetId: 'gke_rect',
});
//
//
// also need to change the name of the output file and folder at end of script

/// this script collects spatial covariates for each transect segment
/// you must import the correct ecosystem extent (rectangle)
/// you must also change the file name at the end to which the covariates are exported 
// the script merges multiple rasters into an image collection

// Modis 500m monthly burned area product from 2001 to 2024
// clipped to viewsheds
var fire = ee.ImageCollection('MODIS/061/MCD64A1')
  .map(function(image){return image.unmask().clip(ecosystem)})
  .filterDate('2001-01-01', '2024-01-01');
var firedate = fire.select(['BurnDate','QA']);

// This function sets the starting position for extracting the quality control bit and the number of
// positions necessary.
var GetQABits = function(image, bitStart, numBits){
    var bitmapped = image.select('QA').rightShift(bitStart).mod(Math.pow(2, numBits));
    return bitmapped;
};
//This function isolates the 2nd bit in Modis QA bitmask which codes for usable vs unusable data
var getUsable = function(image){
  var detailqa = image.select('QA');
  var Qual = GetQABits(detailqa, 1, 1);
  return Qual;
};
// map function across Modis image collection to generate new collection that is only QA bit for usable/unusable data
var Usable = firedate.map(getUsable);
// function to return 1 if burndate is a number (a fire was observed)	
var getfiredays = function(image){
  return image.gte(1);
};
// funciton to return 1 if QA bit indicated usable data
var getnotnulldays = function(image){
  return image.eq(1);
};
// map binary functions over the image collection
var burningdays = firedate.select('BurnDate').map(getfiredays);
var nonnulldays = Usable.map(getnotnulldays);
// for each pixel, sum across months to get total number of months a fire was observed
var sumfiredays = burningdays.sum();//each pixel, how many burn days observed 2000 - 2024

// for each pixel, sum across months to get total number of months with usable data
var sumnotnull = nonnulldays.sum();

// for each pixel, divide number of fires by number of years of clean data
// this translates into burn freqency - burns per year by pixel
// to get years sumtonull is total months of usable data so must be divided by 12 to get fully observed years of data
var freqfire = sumfiredays.divide(sumnotnull.divide(12)).double();
//Map.addLayer(freqfire);


// Land cover
// Load input imagery: Copernicus proportional land cover at 100 m resolution.
var treecover = ee.Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019').select('tree-coverfraction').unmask().double().clip(ecosystem);
var shrubcover = ee.Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019').select('shrub-coverfraction').unmask().double().clip(ecosystem);
var herbaceouscover = ee.Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019').select('grass-coverfraction').unmask().double().clip(ecosystem);
var builtupcover = ee.Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019').select('urban-coverfraction').unmask().double().clip(ecosystem);
var seasonalwatercover = ee.Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019').select('water-seasonal-coverfraction').unmask().double().clip(ecosystem);

//Map.addLayer(treecover);


//Distance to river
var rivers = ee.FeatureCollection("WWF/HydroSHEDS/v1/FreeFlowingRivers");
// Reduce to Zambia - major rivers e.g. Luangwa, Lufupa, etc
var major = rivers.filter('BB_DIS_ORD <= 5').filter('COUNTRY == "Zambia"');
var streams = rivers.filter('BB_DIS_ORD == 6').filter('COUNTRY == "Zambia"');
Map.addLayer(major,{color: 'blue'});
Map.addLayer(streams,{color: 'red'});
// raster of distances to nearest major river
var riverprox =  ee.Image(major.distance(80000)).double().clip(ecosystem);
var streamprox = ee.Image(streams.distance(80000)).double().clip(ecosystem);


// park raster: 1 is park 0 is not park
var parkmask = ee.Image.constant(1);
var park_bnd = parkmask.clip(parks_vec).mask().double();

//gma raster: 1 is gma, 0 is not gma
var gmamask = ee.Image.constant(1);
var gma_bnd = parkmask.clip(gmas_vec).mask().double();
      
//stack all the bands together
var gke_area_covariates = freqfire.addBands(treecover)
.addBands(shrubcover)
.addBands(herbaceouscover)
.addBands(builtupcover)
.addBands(seasonalwatercover)
.addBands(riverprox)
.addBands(streamprox)
.addBands(park_bnd)
.addBands(gma_bnd).rename(['annburnfreq','treecover','shrubcover','herbcover','builtupcover','seasonalwatercover','distriver','diststream','park','gma']);

print(gke_area_covariates);

//Map.addLayer(gke_area_covariates);
Map.addLayer(gke_area_covariates,
    {
        bands: 'treecover',
        palette: ['808000', 'FF0000'],
    },
    'Burn frequency');
    



    
Export.image.toDrive({
  image: gke_area_covariates,
  description: 'lve_area_covariates',
  folder: 'LVE_gee_outputs',
  region: ecosystem,
  scale: 100,
  crs: 'EPSG:4326'
});

