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

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import org.mskcc.cbio.portal.dao.DaoException;
import org.mskcc.cbio.portal.dao.DaoPdbPtmData;
import org.mskcc.cbio.portal.dao.DaoPdbUniprotResidueMapping;
import org.mskcc.cbio.portal.model.PdbUniprotAlignment;
import org.mskcc.cbio.portal.model.PdbUniprotResidueMapping;

/**
 *
 * @author jgao
 */
public class PTM3DHotspotDetective extends ProteinStructureHotspotDetective {

    public PTM3DHotspotDetective(HotspotDetectiveParameters parameters) throws HotspotException {
        super(parameters);
    }
    
    /**
     *
     * @param hotspots
     * @return
     */
    @Override
    protected Map<MutatedProtein,Set<Hotspot>> processSingleHotspotsOnAProtein(MutatedProtein protein,
            Map<Integer, Hotspot> mapResidueHotspot) throws HotspotException {
        if (protein.getProteinLength()>5000) {
            System.out.println("Protein longer than 5000, skipping..");
            return Collections.emptyMap();
        }
        Map<SortedSet<Integer>,Set<Hotspot>> mapResiduesHotspots3D = new HashMap<SortedSet<Integer>,Set<Hotspot>>();
        Map<MutatedProtein3D,List<PdbUniprotAlignment>> mapAlignments = getPdbUniprotAlignments(protein);
        Map<MutatedProtein3D,int[]> alignedRanges = getAlignedRanges(mapAlignments);
        Map<MutatedProtein3D,Map<SortedSet<Integer>,String>> pdbPtms = get3DPTMs(mapAlignments, mapResidueHotspot.keySet());
        
        int largestPos = getLargestAlgnedPosition(alignedRanges);
        if (largestPos > protein.getProteinLength()) {
            protein.setProteinLength(largestPos);
        }
        
        int[] counts = getMutationCountsOnProtein(mapResidueHotspot, protein.getProteinLength());
        
        for (Map.Entry<MutatedProtein3D,Map<SortedSet<Integer>,String>> entryPdbPtm : pdbPtms.entrySet()) {
            MutatedProtein3D protein3D = entryPdbPtm.getKey();
            
            int[][] decoyCountsList = null;
            if (parameters.calculatePValue()) {
                int[] range = alignedRanges.get(protein3D);
                decoyCountsList = generateDecoys(counts, range[0], range[1]+1, 10000);
            }
            
            Map<SortedSet<Integer>,String> ptmMap = entryPdbPtm.getValue();
            for (Map.Entry<SortedSet<Integer>,String> entryPtm : ptmMap.entrySet()) {
                SortedSet<Integer> residues = entryPtm.getKey();
                String ptmLabel = entryPtm.getValue();
                HotspotImpl hotspot3D = new HotspotImpl(protein3D, numberOfsequencedCases, residues);
                hotspot3D.setLabel(ptmLabel);
                
                int maxCap = 0;
                for (int residue : residues) {
                    Hotspot hotspot = mapResidueHotspot.get(residue);
                    hotspot3D.mergeHotspot(hotspot);
                    
                    int num = hotspot.getPatients().size();
                    if (num > maxCap) {
                        maxCap = num;
                    }
                }
                
                if (hotspot3D.getPatients().size()>=parameters.getThresholdSamples()) {
                    Set<Hotspot> hotspots3D = mapResiduesHotspots3D.get(residues);
                    if (hotspots3D==null) {
                        hotspots3D = new HashSet<Hotspot>();
                        mapResiduesHotspots3D.put(residues, hotspots3D);
                    }
                    
                    if (parameters.calculatePValue()) {
                        DetectedInDecoy detectedInDecoy = new StructurePTMHotspotDetectedInDecoy(residues, maxCap, hotspot3D.getPatients().size());
                        double p = getP(detectedInDecoy, decoyCountsList);
                        hotspot3D.setPValue(p);
                    }
                    
                    hotspots3D.add(hotspot3D);
                }
            }
        }
        
        if (mapResiduesHotspots3D.isEmpty()) {
            return Collections.emptyMap();
        }
        
        Set<Hotspot> hotspots3D = new HashSet<Hotspot>(); 
        for (Map.Entry<SortedSet<Integer>,Set<Hotspot>> entryMapResiduesHotspots3D : mapResiduesHotspots3D.entrySet()) {
            SortedSet<Integer> residues = entryMapResiduesHotspots3D.getKey();
            Set<Hotspot> hotspots = entryMapResiduesHotspots3D.getValue();
            Hotspot3D hotspot3D = new Hotspot3DImpl(protein, numberOfsequencedCases, residues, hotspots);
            Hotspot hs = hotspots.iterator().next();
            hotspot3D.mergeHotspot(hs); // add mutations
            hotspot3D.setLabel(hs.getLabel()+" "+hotspot3D.getLabel());
            hotspots3D.add(hotspot3D);
        }
        
        return Collections.singletonMap(protein, hotspots3D);
    }
    
    private int getLargestAlgnedPosition(Map<MutatedProtein3D,int[]> alignedRanges) {
        int pos = Integer.MIN_VALUE;
        for (int[] range : alignedRanges.values()) {
            if (range[1] > pos) {
                pos = range[1];
            }
        }
        
        return pos;
    }
    
    private Map<MutatedProtein3D,int[]> getAlignedRanges(Map<MutatedProtein3D,List<PdbUniprotAlignment>> mapAlignments) {
        Map<MutatedProtein3D,int[]> ranges = new HashMap<MutatedProtein3D,int[]>();
        for (Map.Entry<MutatedProtein3D,List<PdbUniprotAlignment>> entry : mapAlignments.entrySet()) {
            MutatedProtein3D protein3D = entry.getKey();
            List<PdbUniprotAlignment> alignments = entry.getValue();
            int[] range = new int[] {Integer.MAX_VALUE, Integer.MIN_VALUE};
            for (PdbUniprotAlignment alg : alignments) {
                int left = alg.getUniprotFrom();
                int right = alg.getUniprotTo();
                if (left < range[0]) {
                    range[0] = left;
                }
                if (right > range[1]) {
                    range[1] = right;
                }
            }
            ranges.put(protein3D, range);
        }
        return ranges;
    }
    
    private static class StructurePTMHotspotDetectedInDecoy implements DetectedInDecoy {
        private final SortedSet<Integer> residues;
        private final int maxCap;
        private final int targetCount;
        StructurePTMHotspotDetectedInDecoy(final SortedSet<Integer> residues, final int maxCap, final int targetCount) {
            this.residues = residues;
            this.maxCap = maxCap;
            this.targetCount = targetCount;
        }
        public boolean isDetectedInDecoy(final int[] decoy) {
            int count = 0;
            for (int r : residues) {
                int c = decoy[r];
                if (c > maxCap) {
                    c = maxCap;
                }
                count += c;
                if (count >= targetCount) {
                    return true;
                }
            }

            return false;
        }
    }
    
    private Map<MutatedProtein3D,Map<SortedSet<Integer>, String>> get3DPTMs(Map<MutatedProtein3D,List<PdbUniprotAlignment>> mapAlignments, Set<Integer> residues) throws HotspotException {
        try {
            Map<MutatedProtein3D,Map<SortedSet<Integer>, String>> ptms = new HashMap<MutatedProtein3D,Map<SortedSet<Integer>,String>>();
            for (Map.Entry<MutatedProtein3D,List<PdbUniprotAlignment>> entryMapAlignments : mapAlignments.entrySet()) {
                MutatedProtein3D protein3D = entryMapAlignments.getKey();
                List<PdbUniprotAlignment> alignments = entryMapAlignments.getValue();
                
                OneToOneMap<Integer, Integer> pdbUniprotResidueMapping = getPdbUniprotResidueMapping(alignments);
//                protein3D.setProteinLength(pdbUniprotResidueMapping.size()); // only mapped residues
                
                // only retain the mapped mutated residues
                pdbUniprotResidueMapping.retainByValue(residues);
                
                if (pdbUniprotResidueMapping.size()==0) {
                    continue;
                }
                
                Map<SortedSet<Integer>,String> pdbPtms = getPdbPtm(protein3D, pdbUniprotResidueMapping);
                
                ptms.put(protein3D, pdbPtms);
            }
            
            return ptms;
        } catch (DaoException e) {
            throw new HotspotException(e);
        }
    }
    
    protected Map<SortedSet<Integer>, String> getPdbPtm(MutatedProtein3D protein3D, OneToOneMap<Integer, Integer> pdbUniprotResidueMapping) throws DaoException {
        Set<Integer> pdbResiduesAllMapped = pdbUniprotResidueMapping.getKeySet();
        Map<SortedSet<Integer>, String> pdbMap = DaoPdbPtmData.getPdbPtmModules(protein3D.getPdbId(), protein3D.getPdbChain(), pdbResiduesAllMapped);
        Map<SortedSet<Integer>, String> uniprotMap = new HashMap<SortedSet<Integer>,String>();
        for (Map.Entry<SortedSet<Integer>, String> entry : pdbMap.entrySet()) {
            SortedSet<Integer> pdbResidues = entry.getKey();
            String ptm = entry.getValue();
            SortedSet<Integer> uniprotResidues = new TreeSet<Integer>();
            for (Integer pr : pdbResidues) {
                Integer ur = pdbUniprotResidueMapping.getByKey(pr);
                if (ur != null) { // it is possible that some of the pdb residues are not aligned
                    uniprotResidues.add(ur);
                }
            }
            uniprotMap.put(uniprotResidues, ptm);
        }
        return uniprotMap;
    }
    
    private OneToOneMap<Integer, Integer> getPdbUniprotResidueMapping(List<PdbUniprotAlignment> alignments) throws HotspotException {
        Collections.sort(alignments, new Comparator<PdbUniprotAlignment>() {
            @Override
            public int compare(PdbUniprotAlignment align1, PdbUniprotAlignment align2) {
                int ret = align1.getEValue().compareTo(align2.getEValue()); // sort from small evalue to large evalue
                if (ret == 0) {
                    ret = -align1.getIdentityPerc().compareTo(align2.getIdentityPerc());
                }
                return ret;
            }
        });
        
        try {
            OneToOneMap<Integer, Integer> map = new OneToOneMap<Integer, Integer>();
            for (PdbUniprotAlignment alignment : alignments) {
                List<PdbUniprotResidueMapping> mappings = 
                        DaoPdbUniprotResidueMapping.getResidueMappings(alignment.getAlignmentId());
                for (PdbUniprotResidueMapping mapping : mappings) {
                    if (parameters.getIncludingMismatchesFor3DHotspots() || mapping.getMatch().matches("[A-Z]")) { // exact match
                        map.put(mapping.getPdbPos(), mapping.getUniprotPos());
                    }
                }
            }
            return map;
        } catch (DaoException e) {
            throw new HotspotException(e);
        }
    }
}
