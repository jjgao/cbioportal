/** Copyright (c) 2012 Memorial Sloan-Kettering Cancer Center.
**
** This library is free software; you can redistribute it and/or modify it
** under the terms of the GNU Lesser General Public License as published
** by the Free Software Foundation; either version 2.1 of the License, or
** any later version.
**
** This library is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
** MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
** documentation provided hereunder is on an "as is" basis, and
** Memorial Sloan-Kettering Cancer Center 
** has no obligations to provide maintenance, support,
** updates, enhancements or modifications.  In no event shall
** Memorial Sloan-Kettering Cancer Center
** be liable to any party for direct, indirect, special,
** incidental or consequential damages, including lost profits, arising
** out of the use of this software and its documentation, even if
** Memorial Sloan-Kettering Cancer Center 
** has been advised of the possibility of such damage.  See
** the GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this library; if not, write to the Free Software Foundation,
** Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
**/
package org.mskcc.cbio.portal.hotspots;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.codehaus.jackson.map.ObjectMapper;
import org.mskcc.cbio.portal.dao.DaoCancerStudy;
import org.mskcc.cbio.portal.dao.DaoException;
import org.mskcc.cbio.portal.dao.DaoGeneOptimized;
import org.mskcc.cbio.portal.dao.DaoGeneticProfile;
import org.mskcc.cbio.portal.dao.DaoPatient;
import org.mskcc.cbio.portal.dao.DaoSample;
import org.mskcc.cbio.portal.model.CancerStudy;
import org.mskcc.cbio.portal.model.CanonicalGene;
import org.mskcc.cbio.portal.model.ExtendedMutation;
import org.mskcc.cbio.portal.servlet.QueryBuilder;

/**
 *
 * @author jgao
 */
public final class HotspotMain {
    private static Logger logger = Logger.getLogger(HotspotMain.class);
    public static final String HOTSPOT_TYPE = "hotspot_type";
    public static final String MUTATION_TYPE = "mutation_type";
    public static final String PTM_TYPE = "ptm_type";
    public static final String GENES = "genes";
    public static final String THRESHOLD_SAMPLES = "threshold_samples";
    public static final String THRESHOLD_MUTATIONS_HYPERMUTATOR = "threshold_hypermutator";
    public static final String PTM_HOTSPOT_WINDOW = "window_ptm";
    public static final String THRESHOLD_DISTANCE_CLOSEST_DISTANCE_CONTACT_MAP = "threshold_distance_closest_atom_3d";
    public static final String THRESHOLD_DISTANCE_C_ALPHA_CONTACT_MAP = "threshold_distance_c_alpha_3d";
    public static final String THRESHOLD_DISTANCE_ERROR_CONTACT_MAP = "threshold_distance_error_3d";
    public static final String THRESHOLD_UNIPROT_PDB_ALIGNMENT_IDENTP = "threshold_identp_3d";
    public static final String LINEAR_HOTSPOT_WINDOW = "window_linear";
    public static final String SINGLE_HOTSPOT_SEPARATE_BY_PROTEIN_CHANGE = "separate_by_protein_change";
    public static final String THRESHOLD_PREFILTER_RECURRENCE = "threshold_prefilter_recurrence";
    public static final String THRESHOLD_PVALUE = "threshold_pvalue";
    
    public static void main(String[] args) throws IOException, HotspotException  {
        String parameterConfigFile = args[0];
        Map<String, String> map = getRequestParameterMap(parameterConfigFile);
        File file;
        if (args.length>1) {
            file = new File(args[1]);
        } else {
            file = File.createTempFile("hotspots-",".txt");
        }
        System.out.println("Save to: "+file.getAbsolutePath());
        PrintWriter out = new PrintWriter(file);
        try {
            detectHotspot(map, out);
        } finally {
            out.close();
        }
    }
    
    private static Map<String, String> getRequestParameterMap(String parameterConfigFile) throws IOException {
        Properties properties = new Properties();
        properties.load(new FileInputStream(parameterConfigFile));
        
        Map<String, String> map = new HashMap<String, String>();
        for (String key : properties.stringPropertyNames()) {
            map.put(key, properties.getProperty(key));
        }
        return map;
    }
    
    public static void detectHotspot(Map<String,String> requestMap, PrintWriter out)
            throws HotspotException, IOException {
        String studyStableIdsStr = requestMap.get(QueryBuilder.CANCER_STUDY_ID);
        String hotspotType = requestMap.get(HOTSPOT_TYPE);
        String mutationType = requestMap.get(MUTATION_TYPE);
        int threshold = Integer.parseInt(requestMap.get(THRESHOLD_SAMPLES));
        int thresholdHyper = Integer.parseInt(requestMap.get(THRESHOLD_MUTATIONS_HYPERMUTATOR));
        int thresholdPrefilterRecurrence = Integer.parseInt(requestMap.get(THRESHOLD_PREFILTER_RECURRENCE));
        String strThresholdPValue = requestMap.get(THRESHOLD_PVALUE);
        double thresholdPValue = strThresholdPValue==null || strThresholdPValue.isEmpty()?1:Double.parseDouble(strThresholdPValue);
        String genes = requestMap.get(GENES);
        Set<Long>  entrezGeneIds = new HashSet<Long>();
        Set<Long>  excludeEntrezGeneIds = new HashSet<Long>();
        if (genes!=null) {
            DaoGeneOptimized daoGeneOptimized = DaoGeneOptimized.getInstance();
            for (String gene : genes.split("[, ]+")) {
                CanonicalGene canonicalGene = daoGeneOptimized.getGene(gene);
                if (canonicalGene!=null) {
                    entrezGeneIds.add(canonicalGene.getEntrezGeneId());
                } else if (gene.startsWith("-")) {
                    canonicalGene = daoGeneOptimized.getGene(gene.substring(1));
                    if (canonicalGene!=null) {
                        excludeEntrezGeneIds.add(canonicalGene.getEntrezGeneId());
                    }
                }
            };
        }
        
        Set<Hotspot> hotspots = Collections.emptySet();
        Map<Integer,String> cancerStudyIdMapping = new HashMap<Integer,String>();
        String[] studyStableIds = studyStableIdsStr.split("[, ]+");
        
        Set<Integer> studyIds = new HashSet<Integer>();
        for (String stableId : studyStableIds) {
            CancerStudy study = null;
            try {
                study = DaoCancerStudy.getCancerStudyByStableId(stableId);
            } catch (DaoException e) {
                throw new HotspotException(e);
            }
            if (study!=null) {
                studyIds.add(study.getInternalId());
                cancerStudyIdMapping.put(study.getInternalId(), stableId);
            }
        }

        HotspotDetectiveParameters hotspotDetectiveParameters = new HotspotDetectiveParametersImpl();
        hotspotDetectiveParameters.setCancerStudyIds(studyIds);
        hotspotDetectiveParameters.setEntrezGeneIds(entrezGeneIds);
        hotspotDetectiveParameters.setExcludeEntrezGeneIds(excludeEntrezGeneIds);
        hotspotDetectiveParameters.setMutationTypes(Arrays.asList(mutationType.split("[, ]+")));
        hotspotDetectiveParameters.setThresholdHyperMutator(thresholdHyper);
        hotspotDetectiveParameters.setThresholdSamples(threshold);
        hotspotDetectiveParameters.setPrefilterThresholdSamplesOnSingleResidue(thresholdPrefilterRecurrence);
        hotspotDetectiveParameters.setPValueThreshold(thresholdPValue);

        HotspotDetective hotspotDetective;
        if (hotspotType.equalsIgnoreCase("single")) {
            boolean separateByProteinChange = Boolean.parseBoolean(requestMap.get(SINGLE_HOTSPOT_SEPARATE_BY_PROTEIN_CHANGE));
            hotspotDetectiveParameters.setSeperateByProteinChangesForSingleResidueHotspot(separateByProteinChange);
            hotspotDetective = new SingleHotspotDetective(hotspotDetectiveParameters);
        } else if (hotspotType.equalsIgnoreCase("linear")) {
            int window = Integer.parseInt(requestMap.get(LINEAR_HOTSPOT_WINDOW));
            hotspotDetectiveParameters.setLinearSpotWindowSize(window);
            hotspotDetective = new LinearHotspotDetective(hotspotDetectiveParameters);
        } else if (hotspotType.equalsIgnoreCase("3d")) {
            String strThresholdDisClosestAtoms = requestMap.get(THRESHOLD_DISTANCE_CLOSEST_DISTANCE_CONTACT_MAP);
            double thresholdDisClosestAtoms = strThresholdDisClosestAtoms==null?0:Double.parseDouble(strThresholdDisClosestAtoms);
            hotspotDetectiveParameters.setDistanceClosestAtomsThresholdFor3DHotspots(thresholdDisClosestAtoms);

            String strThresholdDisCAlpha = requestMap.get(THRESHOLD_DISTANCE_C_ALPHA_CONTACT_MAP);
            double thresholdDisCAlpha = strThresholdDisCAlpha==null?0:Double.parseDouble(strThresholdDisCAlpha);
            hotspotDetectiveParameters.setDistanceCAlphaThresholdFor3DHotspots(thresholdDisCAlpha);

            String strThresholdDisError = requestMap.get(THRESHOLD_DISTANCE_ERROR_CONTACT_MAP);
            double thresholdDisError = strThresholdDisError==null?0:Double.parseDouble(strThresholdDisError);
            hotspotDetectiveParameters.setDistanceErrorThresholdFor3DHotspots(thresholdDisError);



            String strThresholdIdentp = requestMap.get(THRESHOLD_UNIPROT_PDB_ALIGNMENT_IDENTP);
            double thresholdIdentp = strThresholdIdentp==null?0:Double.parseDouble(strThresholdIdentp);
            hotspotDetectiveParameters.setIdentpThresholdFor3DHotspots(thresholdIdentp);
            hotspotDetectiveParameters.setIncludingMismatchesFor3DHotspots(false);

            hotspotDetective = new ProteinStructureHotspotDetective(hotspotDetectiveParameters);
        } else if (hotspotType.startsWith("ptm")) {
            int thresholdDis = Integer.parseInt(requestMap.get(PTM_HOTSPOT_WINDOW));
//                String ptmType = requestMap.get(PTM_TYPE);
            hotspotDetectiveParameters.setPtmHotspotWindowSize(thresholdDis);
            hotspotDetective = new PTMHotspotDetective(hotspotDetectiveParameters);
//            } else if (type.equalsIgnoreCase("truncating-sep")) {
//                 mapKeywordStudyCaseMut = DaoMutation.getTruncatingMutatationStatistics(
//                    studyIds.toString(), threshold, concatEntrezGeneIds, concatExcludeEntrezGeneIds);
        } else if (hotspotType.equalsIgnoreCase("3d-ptm")) {
            String strThresholdIdentp = requestMap.get(THRESHOLD_UNIPROT_PDB_ALIGNMENT_IDENTP);
            double thresholdIdentp = strThresholdIdentp==null?0:Double.parseDouble(strThresholdIdentp);
            hotspotDetectiveParameters.setIdentpThresholdFor3DHotspots(thresholdIdentp);
            hotspotDetectiveParameters.setIncludingMismatchesFor3DHotspots(false);

             hotspotDetective = new PTM3DHotspotDetective(hotspotDetectiveParameters);
        } else {
            throw new IllegalStateException("wrong hotspot type: "+hotspotType);
        }

        hotspotDetective.detectHotspot();
        hotspots = hotspotDetective.getDetectedHotspots();
        
        // transform the data to use stable cancer study id
        Map<String,Map<String, Map<String,Set<String>>>> map =
                new HashMap<String,Map<String, Map<String,Set<String>>>>(hotspots.size());
        for (Hotspot hotspot : hotspots) {
            String label = hotspot.getLabel();
            Map<String, Map<String,Set<String>>> map1 = new HashMap<String, Map<String,Set<String>>>();
            for (ExtendedMutation mutation : hotspot.getMutations()) {
                String cancerStudy = null;
                try {
                    cancerStudy = DaoCancerStudy.getCancerStudyByInternalId(
                        DaoGeneticProfile.getGeneticProfileById(
                        mutation.getGeneticProfileId()).getCancerStudyId()).getCancerStudyStableId();
                } catch (DaoException e) {
                    throw new HotspotException(e);
                }
                Map<String,Set<String>> map2 = map1.get(cancerStudy);
                if (map2==null) {
                    map2 = new HashMap<String,Set<String>>();
                    map1.put(cancerStudy, map2);
                }
                
                int sampleId = mutation.getSampleId();
                String patientId = DaoPatient.getPatientById(DaoSample.getSampleById(sampleId).getInternalPatientId()).getStableId();
                
                Set<String> aaChanges = map2.get(patientId);
                if (aaChanges==null) {
                    aaChanges = new HashSet<String>();
                    map2.put(patientId, aaChanges);
                }
                aaChanges.add(mutation.getProteinChange());
            }
            map.put(label, map1);
        }

        String format = requestMap.get("format");
        
        if (format==null || format.equalsIgnoreCase("json")) {
            ObjectMapper mapper = new ObjectMapper();
            out.write(mapper.writeValueAsString(map));
        } else if (format.equalsIgnoreCase("text")) {
            out.write("Alteration\t");
            out.write(StringUtils.join(studyStableIds,"\t"));
            out.write("\n");
            for (Map.Entry<String,Map<String, Map<String,Set<String>>>> entry : map.entrySet()) {
                String keyword = entry.getKey();
                out.write(keyword);
                Map<String, Map<String,Set<String>>> mapStudyCaseMut = entry.getValue();
                for (String study : studyStableIds) {
                    Map<String,Set<String>> mapCaseMut = mapStudyCaseMut.get(study);
                    out.write("\t");
                    if (mapCaseMut!=null && !mapCaseMut.isEmpty()) {
                        out.write(Integer.toString(mapCaseMut.size()));
                    }
                }
                out.write("\n");
            }
        } 
    }

    private HotspotMain() {
    }
}
