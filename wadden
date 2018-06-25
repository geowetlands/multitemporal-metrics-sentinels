/*////////////////////////////////////////////////////////////////
Creates multitemporal composites out of Sentinel-1 and -2 imagery
Script designed to support wetland mapping by calculating temporal metrics
Author: Javier Muro (javi.estonia@gmail.com)
////////////////////////////////////////////////////////////////*/

/////////////////////////////////Set search parameters///////////////////////////
//Wadden sea
var geom = ee.FeatureCollection('ft:1rhyhRa-UYJJgGtr8PNrW3vywj7Vf2ix0Vzr5cMUh');

var startDate = ee.Date('2016-06-01')
var stopDate = ee.Date('2018-06-01')

////////////////////////////////Auxiliary functions/////////////////////////////
//Visualization parameters
var imageVisParam = {bands: ["NDBI_p90","NDVI_p90","NDWI_p90"],
gamma: 1,
max: 0.7910480499267578,
min: -0.3163292109966278,
opacity: 1
}

var imageVisParam2 = {bands: ["VV_p95","VV_p50","VV_p5"],
gamma: 1,
max: -1.1474982500076294,
min: -56.507240295410156,
opacity: 1
}

var imageVisParam3 = {bands: ["B8_p50","B5_p50","red_p50"],
gamma: 1,
max: 0.36739999055862427,
min: 0.029200000688433647,
opacity: 1
}

//Index Sentinel-2 imagery
function sentinel2toa(img) {
  var toa = img.select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10', 'B11','B12'],  
                       ['aerosol','blue','green','red','B5','red2','B7','B8','red4','h2o','cirrus','swir1','swir2'])
                       .divide(10000)
                       .addBands(img.select(['QA60']))
                       .set('solar_azimuth',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                       .set('solar_zenith',img.get('MEAN_SOLAR_ZENITH_ANGLE'))
    return toa
}

//Function to mask clouds
function ESAcloud(toa) {
  // author: Nick Clinton
  var qa = toa.select('QA60');
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);
  
  // clear if both flags set to zero.
  var clear = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  var cloud = clear.eq(0)
  return cloud
}

function shadowMask(toa,cloud){
  // Author: Gennadii Donchyts
  // License: Apache 2.0
  
  // solar geometry (radians)
  var azimuth =ee.Number(toa.get('solar_azimuth')).multiply(Math.PI).divide(180.0).add(ee.Number(0.5).multiply(Math.PI));
  var zenith  =ee.Number(0.5).multiply(Math.PI ).subtract(ee.Number(toa.get('solar_zenith')).multiply(Math.PI).divide(180.0));

  // find where cloud shadows should be based on solar geometry
  var nominalScale = cloud.projection().nominalScale();
  var cloudHeights = ee.List.sequence(200,10000,500);
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    var shadowVector = zenith.tan().multiply(cloudHeight);
    var x = azimuth.cos().multiply(shadowVector).divide(nominalScale).round();
    var y = azimuth.sin().multiply(shadowVector).divide(nominalScale).round();
    return cloud.changeProj(cloud.projection(), cloud.projection().translate(x, y));
  });
  var potentialShadow = ee.ImageCollection.fromImages(shadows).max();
  
  // shadows are not clouds
  var potentialShadow = potentialShadow.and(cloud.not());
  
  // (modified by Sam Murphy) dark pixel detection 
  var darkPixels = toa.normalizedDifference(['green', 'swir2']).gt(0.25).rename(['dark_pixels']);
  
  // shadows are dark
  var shadow = potentialShadow.and(darkPixels).rename('shadows');
  return shadow
}

function cloud_and_shadow_mask(img) {
  var toa = sentinel2toa(img)
  var cloud = ESAcloud(toa)
  var shadow = shadowMask(toa,cloud)
  var mask = cloud.or(shadow).eq(0)
  return toa.updateMask(mask)
}

//Functions to calculate normalized difference indices
var addNDVI = function(image) {
  return image.addBands(image.normalizedDifference(['red4','red']).rename('NDVI'))
};
var addNDWI = function(image) {
  return image.addBands(image.normalizedDifference(['green','red4']).rename('NDWI'))
};
var addNDBI = function(image) {
  return image.addBands(image.normalizedDifference(['swir1','red4']).rename('NDBI'))
};


//Functions to mask noise at the edges
function maskEdge(img) {
  var mask = img.select(0).unitScale(-25, 10).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
  return img.updateMask(mask.select(0));  
}

//Additional optional function to mask edges. Apply only if needed
var maskedge2 = function (img) {
  return img.clip(img.geometry().buffer(-4000));
};

//////////////////////////////End of Auxiliary functions//////////////////////////

//Load Sentinel-2 collection and apply masks
var images = ee.ImageCollection('COPERNICUS/S2')
  .filterBounds(geom)
  .filterDate(startDate, stopDate)

//Apply cloud mask and add normalized indices
var masked_images = images
  .map(cloud_and_shadow_mask)
  .map(addNDBI)
  .map(addNDVI)
  .map(addNDWI)

//print(masked_images);
  
var NDVIperc90 = masked_images.select('NDVI').reduce(ee.Reducer.percentile({percentiles: [90]}));
var NDWIperc90 = masked_images.select('NDWI').reduce(ee.Reducer.percentile({percentiles: [90]}));
var NDBIperc90 = masked_images.select('NDBI').reduce(ee.Reducer.percentile({percentiles: [90]}));
var NDVIperc50 = masked_images.select('NDVI').reduce(ee.Reducer.percentile({percentiles: [50]}));
var NDWIperc50 = masked_images.select('NDWI').reduce(ee.Reducer.percentile({percentiles: [50]}));
var NDBIperc50 = masked_images.select('NDBI').reduce(ee.Reducer.percentile({percentiles: [50]}));

//print(NDVIperc90)
//print(NDWIperc90)
//print(NDBIperc90)
//print(NDVIperc50)
//print(NDWIperc50)
//print(NDBIperc50)

var B1perc50 = masked_images.select('aerosol').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B2perc50 = masked_images.select('blue').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B3perc50 = masked_images.select('green').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B4perc50 = masked_images.select('red').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B5perc50 = masked_images.select('B5').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B6perc50 = masked_images.select('red2').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B7perc50 = masked_images.select('B7').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B8perc50 = masked_images.select('B8').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B8Aperc50 = masked_images.select('red4').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B9perc50 = masked_images.select('h2o').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B10perc50 = masked_images.select('cirrus').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B11perc50 = masked_images.select('swir1').reduce(ee.Reducer.percentile({percentiles: [50]}));
var B12perc50 = masked_images.select('swir2').reduce(ee.Reducer.percentile({percentiles: [50]}));


// Load the Sentinel-1 ImageCollection
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD');

// Filter by metadata properties.
var vh = sentinel1
  // Filter to get images with VV and VH dual polarization.
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
//  .filter(ee.Filter.eq('platform_number', 'A'))
  .filter(ee.Filter.eq('relativeOrbitNumber_start', 37)) //this correspond to the Wadden sea. Other paths create noise at seams between images
  .filterDate(startDate, stopDate)
  .filterBounds(geom)
;
        


// Filter to get images from different look angles.
var vhAscending = vh
  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
  .map(maskEdge)
  //.map(maskedge2)  //optional
  ;
var vhDescending = vh
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .map(maskEdge)
  //.map(maskedge2)  //optional
  ;

//print(vhDescending)
//print(vhAscending)

//Create a composite at different polarizations. 
//We choose descending because it has less seams
//VH is for some reason wrong in the western swadth over the wadden sea
//Standard deviations have very small seams between images that cannot be removed
var compositeS1 = ee.Image.cat([
  vhDescending.select('VV').reduce(ee.Reducer.percentile({percentiles: [95]})),
  vhDescending.select('VV').reduce(ee.Reducer.percentile({percentiles: [50]})),
  vhDescending.select('VV').reduce(ee.Reducer.percentile({percentiles: [5]})),
  vhDescending.select('VV').reduce(ee.Reducer.stdDev()),
  vhDescending.select('VH').reduce(ee.Reducer.percentile({percentiles: [95]})),
  vhDescending.select('VH').reduce(ee.Reducer.percentile({percentiles: [50]})),
  vhDescending.select('VH').reduce(ee.Reducer.percentile({percentiles: [5]})),
  vhDescending.select('VH').reduce(ee.Reducer.stdDev())
  ]);
print(compositeS1, 'Sentinel-1 composite');

var indicesStack90_50 = ee.Image.cat([NDBIperc90, NDVIperc90, NDWIperc90, NDBIperc50, NDVIperc50, NDWIperc50
                               ]).clip(geom)
print(indicesStack90_50, 'Indices composite')

var sentinel2Stack50 = ee.Image.cat([B2perc50, B3perc50, B4perc50, B5perc50, B6perc50, B7perc50, B8perc50, B8Aperc50, B11perc50, B12perc50,
                               ]).clip(geom)
print(sentinel2Stack50, 'Sentinel-2 median composite')

//Visualize and check results
Map.addLayer(compositeS1.clip(geom), imageVisParam2, 'S1')
Map.addLayer(indicesStack90_50.clip(geom), imageVisParam, 'S2 NDIs')
Map.addLayer(sentinel2Stack50.clip(geom), imageVisParam3, 'S2')



///////////////////////////////////////////////////////////////Exports///////////////////////////////////////////
//Export composites to asset for classification
Export.image.toAsset({
  image: indicesStack90_50.clip(geom).select('NDBI_p90', 'NDVI_p90', 'NDWI_p90', 'NDBI_p50', 'NDVI_p50', 'NDWI_p50').float(), 
  description: 'NDBVWI_90_50_wadden',
  scale: 20,
  crs: 'EPSG:25831',
  region: geom.geometry().bounds(),
  maxPixels: 6000000000,
});

Export.image.toAsset({
  image: sentinel2Stack50.clip(geom).select('blue_p50', 'green_p50', 'red_p50', 'B5_p50', 'red2_p50', 'B7_p50', 'B8_p50', 'red4_p50', 'swir1_p50', 'swir2_p50')
                                    .float(), 
  description: 'S2_perc50_wadden',
  scale: 20,
  crs: 'EPSG:25831',
  region: geom.geometry().bounds(),
  maxPixels: 6000000000,
});

Export.image.toAsset({
  image: compositeS1.clip(geom).select('VV_p95', 'VV_p50', 'VV_p5', 'VH_p95', 'VH_p50', 'VH_p5', 'VV_stdDev', 'VH_stdDev').float(), 
  description: 'S1_VV_VH_95_50_5_SDs_wadden',
  scale: 20,
  crs: 'EPSG:25831',
  region: geom.geometry().bounds(),
  maxPixels: 6000000000,
});
///////////////////////////////////////////////////////////////////
