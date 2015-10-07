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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import org.apache.commons.lang.ArrayUtils;
import org.mskcc.cbio.portal.dao.DaoException;
import org.mskcc.cbio.portal.dao.DaoPdbUniprotResidueMapping;
import org.mskcc.cbio.portal.dao.DaoProteinContactMap;
import org.mskcc.cbio.portal.model.PdbUniprotAlignment;
import org.mskcc.cbio.portal.model.PdbUniprotResidueMapping;

/**
 *
 * @author jgao
 */
public class ProteinStructureHotspotDetective extends AbstractHotspotDetective {

    public ProteinStructureHotspotDetective(HotspotDetectiveParameters parameters) throws HotspotException {
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
        
        List<Integer> counts = getMutationCountsOnProtein(mapResidueHotspot, protein.getProteinLength());
        List<Integer>[] decoyCountsList = generateDecoys(counts, 1000);
        
        Map<SortedSet<Integer>,Set<Hotspot>> mapResiduesHotspots3D = new HashMap<SortedSet<Integer>,Set<Hotspot>>();
        Map<MutatedProtein3D,boolean[][]> contactMaps = getContactMaps(protein);
        int i = 0;
        for (Map.Entry<MutatedProtein3D, boolean[][]> entryContactMaps : contactMaps.entrySet()) {
            MutatedProtein3D protein3D = entryContactMaps.getKey();
            boolean[][] contactMap = entryContactMaps.getValue();
            
            System.out.println("\t"+(++i)+"/"+contactMaps.size()+". Processing "+protein3D.getPdbId()+"."+protein3D.getPdbChain());
            
            Set<SortedSet<Integer>> clusters = findConnectedNeighbors(contactMap, mapResidueHotspot.keySet());
            for (SortedSet<Integer> residues : clusters) {
                if (residues.size()<=1) {
                    continue;
                }
                
                Hotspot hotspot3D = new HotspotImpl(protein3D, numberOfsequencedCases, residues);
                
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
                    
                    double p = getP(contactMap, decoyCountsList, maxCap, hotspot3D.getPatients().size());
                    hotspot3D.setPValue(p);
                    
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
            hotspot3D.mergeHotspot(hotspots.iterator().next()); // add mutations
            hotspots3D.add(hotspot3D);
        }
        
        return Collections.singletonMap(protein, hotspots3D);
    }
    
    private List<Integer> getMutationCountsOnProtein(Map<Integer, Hotspot> mapResidueHotspot, int len) {
        int[] ret = new int[len+1];
        for (Map.Entry<Integer, Hotspot> entry : mapResidueHotspot.entrySet()) {
            ret[entry.getKey()] = entry.getValue().getPatients().size();
        }
        return Arrays.asList(ArrayUtils.toObject(ret));
    }
    
    private List<Integer>[] generateDecoys(List<Integer> counts, int times) {
        List<Integer>[] decoys = new List[times];
        for (int i=0; i<times; i++) {
            decoys[i] = new ArrayList<Integer>(counts);
            Collections.shuffle(decoys[i]);
        }
        return decoys;
    }
    
    private double getP(boolean[][] graph, List<Integer>[] decoyCountsList, int maxCap, int targetCount) {
        int d = 0;
        for (int i=0; i<decoyCountsList.length; i++) {
            if (isDetectedInDecoy(graph, decoyCountsList[i], maxCap, targetCount)) {
                d++;
            }
        }
        
        return 1.0 * d / decoyCountsList.length;
    }
    
    private boolean isDetectedInDecoy(boolean[][] graph, List<Integer> decoyCounts, int maxCap, int targetCount) {
        for (int i=1; i<graph.length; i++) {
            int count = 0;
            for (int j=1; j<graph.length; j++) {
                if (graph[i][j]) {
                    int c = decoyCounts.get(j);
                    if (c > maxCap) {
                        c = maxCap;
                    }
                    count += c;
                    if (count >= targetCount) {
                        return true;
                    }
                }
            }
        }
        
        return false;
    }
    
    private Set<SortedSet<Integer>> findConnectedNeighbors(boolean[][] graph, Set<Integer> nodes) {
        Set<SortedSet<Integer>> ret = new HashSet<SortedSet<Integer>>(nodes.size());
        for (Integer node : nodes) {
            SortedSet<Integer> set = new TreeSet<Integer>();
            set.add(node);
            for (Integer neighbor : nodes) {
                if (graph[node][neighbor]) {
                    set.add(neighbor);
                }
            }
            ret.add(set);
        }
        return ret;
    }
    
    @Override
    protected int getLengthOfProtein(MutatedProtein protein, Collection<Hotspot> mapResidueHotspot) {
        return super.getLengthOfProtein(protein, mapResidueHotspot);
        //return protein.getProteinLength(); // skip it.. since it's already set
    }
    
    private List<SortedSet<Integer>> cleanResiduesSet(Collection<SortedSet<Integer>> residuesSet) {
        List<SortedSet<Integer>> ret = new ArrayList<SortedSet<Integer>>();
        for (SortedSet<Integer> residues : residuesSet) {
            boolean exist = false;
            for (SortedSet<Integer> existing : ret) {
                if (existing.containsAll(residues)) {
                    exist = true;
                    break;
                }
                
                if (residues.containsAll(existing)) {
                    existing.addAll(residues);
                    exist = true;
                    break;
                }
            }
            
            if (!exist) {
                ret.add(residues);
            }
        }
        return ret;
    }
    
    private Map<MutatedProtein3D,boolean[][]> getContactMaps(MutatedProtein protein) throws HotspotException {
        try {
            Map<MutatedProtein3D,boolean[][]> contactMaps = new HashMap<MutatedProtein3D,boolean[][]>();
            Map<MutatedProtein3D,List<PdbUniprotAlignment>> mapAlignments = getPdbUniprotAlignments(protein);
            for (Map.Entry<MutatedProtein3D,List<PdbUniprotAlignment>> entryMapAlignments : mapAlignments.entrySet()) {
                MutatedProtein3D protein3D = entryMapAlignments.getKey();
                List<PdbUniprotAlignment> alignments = entryMapAlignments.getValue();
                OneToOneMap<Integer, Integer> pdbUniprotResidueMapping = getPdbUniprotResidueMapping(alignments);
                Map<Integer, Set<Integer>> pdbContactMap = getPdbContactMap(protein3D);
                
                boolean[][] contactMap = getUniProtContactMap(pdbContactMap, pdbUniprotResidueMapping, protein.getProteinLength());
                
                contactMaps.put(protein3D, contactMap);
            }
            
            return contactMaps;
        } catch (DaoException e) {
            throw new HotspotException(e);
        }
    }
    
    private boolean[][] getUniProtContactMap(Map<Integer, Set<Integer>> pdbContactMap, OneToOneMap<Integer, Integer> pdbUniprotResidueMapping, int proteinLengh) {
        boolean[][] ret = new boolean[proteinLengh+1][proteinLengh+1];
        
        for (Map.Entry<Integer, Set<Integer>> entryPdbContactMap : pdbContactMap.entrySet()) {
            int pdbPos = entryPdbContactMap.getKey();
            if (!pdbUniprotResidueMapping.hasKey(pdbPos)) {
                continue;
            }
            
            int uniprotPos = pdbUniprotResidueMapping.getByKey(pdbPos);
            if (uniprotPos > proteinLengh) {
                System.err.println("UniProt length longer than protein length");
                continue;
            }
            Set<Integer> pdbNeighbors = entryPdbContactMap.getValue();
            if (pdbNeighbors.isEmpty()) {
                continue;
            }

            for (Integer pdbNeighbor : pdbNeighbors) {
                if (!pdbUniprotResidueMapping.hasKey(pdbNeighbor)) {
                    continue;
                }
                int uniprotNeighbor = pdbUniprotResidueMapping.getByKey(pdbNeighbor);
                if (uniprotNeighbor > proteinLengh) {
                    System.err.println("UniProt length longer than protein length");
                    continue;
                }
                ret[uniprotPos][uniprotNeighbor] = true;
                ret[uniprotNeighbor][uniprotPos] = true;
            }
        }
        
        return ret;
    }
    
    private Map<Integer, Set<Integer>> getPdbContactMap(MutatedProtein3D protein3D) throws DaoException {
        return DaoProteinContactMap.getProteinContactMap(protein3D.getPdbId(), protein3D.getPdbChain(),
                        null, parameters.getDistanceClosestAtomsThresholdFor3DHotspots(),
                        parameters.getDistanceCAlphaThresholdFor3DHotspots(),
                        parameters.getDistanceErrorThresholdFor3DHotspots());
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
    
    protected Map<MutatedProtein3D,List<PdbUniprotAlignment>> getPdbUniprotAlignments(MutatedProtein protein) throws HotspotException {
        try {
            List<PdbUniprotAlignment> alignments = DaoPdbUniprotResidueMapping.getAlignments(protein.getUniprotId(),
                    Double.valueOf(parameters.getIdentpThresholdFor3DHotspots()).floatValue());
            Map<MutatedProtein3D,List<PdbUniprotAlignment>> map = new HashMap<MutatedProtein3D,List<PdbUniprotAlignment>>();
            for (PdbUniprotAlignment alignment : alignments) {
                MutatedProtein3DImpl protein3D = new MutatedProtein3DImpl(protein);
                protein3D.setPdbId(alignment.getPdbId());
                protein3D.setPdbChain(alignment.getChain());
                
                List<PdbUniprotAlignment> list = map.get(protein3D);
                if (list==null) {
                    list = new ArrayList<PdbUniprotAlignment>();
                    map.put(protein3D, list);
                }
                list.add(alignment);
            }
            return map;
        } catch (DaoException e) {
            throw new HotspotException(e);
        }
    }
    
    protected final class OneToOneMap<K, V> {
        private Map<K, V> keyToValMap;
        private Map<V, K> valToKeyMap;
        
        OneToOneMap() {
            keyToValMap = new HashMap<K,V>();
            valToKeyMap = new HashMap<V,K>();
        }
        
        void put(K k, V v) {
            if (!keyToValMap.containsKey(k) && !valToKeyMap.containsKey(v)) {
                keyToValMap.put(k, v);
                valToKeyMap.put(v, k);
            }
        }
        
        int size() {
            return keyToValMap.size();
        }
        
        boolean hasKey(K k) {
            return keyToValMap.containsKey(k);
        }
        
        boolean hasValue(V v) {
            return valToKeyMap.containsKey(v);
        }
        
        V getByKey(K k) {
            return keyToValMap.get(k);
        }
        
        K getByValue(V v) {
            return valToKeyMap.get(v);
        }

        Set<K> getKeySet() {
            return keyToValMap.keySet();
        }

        Set<V> getValSet() {
            return valToKeyMap.keySet();
        }
        
        void retainByKey(Collection<K> keys) {
            keyToValMap.keySet().retainAll(keys);
            Iterator<Map.Entry<V,K>> itEntry = valToKeyMap.entrySet().iterator();
            while (itEntry.hasNext()) {
                Map.Entry<V,K> entry = itEntry.next();
                if (!keyToValMap.containsKey(entry.getValue())) {
                    itEntry.remove();
                }
            }
        }
        
        void retainByValue(Collection<V> values) {
            valToKeyMap.keySet().retainAll(values);
            Iterator<Map.Entry<K,V>> itEntry = keyToValMap.entrySet().iterator();
            while (itEntry.hasNext()) {
                Map.Entry<K,V> entry = itEntry.next();
                if (!valToKeyMap.containsKey(entry.getValue())) {
                    itEntry.remove();
                }
            }
        }
    }
}
