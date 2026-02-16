// ================================================================
// Section 2.3 Preparation of Sentinel-2-derived vegetation indices
// ================================================================

/* The computation of spectral indices will be followed by Montero et al. (2022) 
from https://github.com/awesome-spectral-indices/spectral */

// STEP 1: REQUIRE THE SPECTRAL MODULE
var spectral = require("users/dmlmont/spectral:spectral");

// STEP 2: DEFINE THE STUDY SITE
/* The study site was at Nongbua Sa-Art Sub-distrinc, Buayai District, Nakhon Ratchasrima Province.
There were 47 cassava fields on the study site */

var nongbua = ee.Geometry.Polygon(
        [[[102.24965, 15.58734],
          [102.24965, 15.5596],
          [102.27575, 15.5596],
          [102.27575, 15.58734]]], null, false);
          
Map.addLayer(nongbua, {color:'white'}, 'Nongbua');
Map.centerObject(nongbua, 14);
Map.setOptions('SATELLITE');

/* There were 47 cassava fields of the study site */
var nongbua47 = ee.FeatureCollection('users/kanokpornpromn/Nongbua_47plots');
Map.addLayer(nongbua47, {color:'white'}, 'Nongbua47plots');
Map.centerObject(nongbua47, 14);
Map.setOptions('SATELLITE');

// STEP 3: DEFINE THE DATES TO FILTER - THE DATE THAT HAVE NO-CLOUD
var datesToFilter = [
  '2022-01-11',
  '2022-01-16',
  '2022-01-26',
  '2022-02-25',
  '2022-03-12',
  '2022-03-17',
  '2022-03-27',
  '2022-04-11',
  '2022-04-21',
  '2022-04-26',
  '2022-06-20'
];

// STEP 4: REQUIRE THE DATASET: SENTINEL-2 SR HARMONIZED Level-2A
var dataset = 'COPERNICUS/S2_SR_HARMONIZED';

// Create an empty image collection to store the filtered images
var filteredCollection = ee.ImageCollection([]);

// Function to mask clouds using the QA60 band
function maskS2clouds(img) {
  var qa = img.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return img.updateMask(mask).divide(10000);
}

//Function to map over an image colletion
function addIndices(img) {
  
  // First apply cloud mask
  img = maskS2clouds(img);
  
  // Required parameters according to the required bands
  var parameters = {
    "N": img.select("B8"),
    "R": img.select("B4"),
    "N2": img.select("B8A"),
    "RE1": img.select("B5"),
    "B": img.select("B2"),
    "G": img.select("B3"),
    "S1": img.select("B11"),
    "S2": img.select("B12"),
    "RE3": img.select("B7"),
    "RE2": img.select("B6"),
    "L": 1.0,
    "g" : 2.5,
    "C1" : 6,
    "C2" : 7.5,
    "C" : 1,
    "A" : 3,
    "D" : 0.2
  };
  
  // scale the image
  img = spectral.scale(img,dataset);

  // compute the new bands of vegetation indices
  return spectral.computeIndex(img,['SeLI','GNDVI','NDWI','NDVI','SAVI','DVI',
  'GRVI','CIG','BNDVI','EVI','RVI','VIG','TCARI'],parameters);
}

// Loop through the dates and filter the collection
datesToFilter.forEach(function(date) {
  var startDate = ee.Date(date);
  var endDate = startDate.advance(1, 'day');
  var filtered = ee.ImageCollection(dataset)
    .select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12',
    'AOT','TCI_R','TCI_G','TCI_B','MSK_CLDPRB','QA10','QA20','QA60'])
    .filterDate(startDate, endDate)
    .filterBounds(nongbua)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(addIndices);
  filteredCollection = filteredCollection.merge(filtered);
});

// STEP 5: PRINT IMAGE COLLECTION HAVING OUR NEW BANDS
print(filteredCollection, 'Sen2_CloudMask having our new bands');

// STEP 6: IMPORT THE BATCH EXPORT FUNCTION AND DOWNLOAD THE IMAGERY WITH 13-VIS BANDS
var batch = require('users/fitoprincipe/geetools:batch');

batch.Download.ImageCollection.toDrive(filteredCollection, 'Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2021-10-01_2023-01-01', 
                {description: 'Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2021-10-01_2023-01-01',
                //filenameprefix: 'Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_',
                scale: 10, // ADJUST TO DESIRED PIXEL RESOLUTION
                region: nongbua,
                fileformat: 'GEOTIFF', 
                crs: 'EPSG:4326',
                type: 'float'});
