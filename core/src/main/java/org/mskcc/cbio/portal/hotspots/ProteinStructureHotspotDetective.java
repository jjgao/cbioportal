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
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
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
        
        int[] counts = getMutationCountsOnProtein(mapResidueHotspot, protein.getProteinLength());
        
        Map<SortedSet<Integer>,Set<Hotspot>> mapResiduesHotspots3D = new HashMap<SortedSet<Integer>,Set<Hotspot>>();
        Map<MutatedProtein3D,ContactMap> contactMaps = getContactMaps(protein);
        int i = 0;
        for (Map.Entry<MutatedProtein3D, ContactMap> entryContactMaps : contactMaps.entrySet()) {
            MutatedProtein3D protein3D = entryContactMaps.getKey();
            ContactMap contactMap = entryContactMaps.getValue();
            
            if (contactMap.getProteinRight()>protein.getProteinLength()) {
                System.err.println("\tMapped Protein resisue longer than protein length.");
                continue;
            }
            
            int[][] decoyCountsList = null;
            
            if (parameters.calculatePValue()) {
                decoyCountsList = generateDecoys(counts, contactMap.getProteinLeft(), contactMap.getProteinRight()+1, 10000);
            }
            
            System.out.println("\t"+protein3D.getGene().getHugoGeneSymbolAllCaps()+" "+(++i)+"/"+contactMaps.size()+". Processing "
                    +protein3D.getPdbId()+"."+protein3D.getPdbChain());
            
            Set<SortedSet<Integer>> clusters = findConnectedNeighbors(contactMap.getContact(), mapResidueHotspot.keySet());
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
                    
                    if (parameters.calculatePValue()) {
                        DetectedInDecoy detectedInDecoy = new StructureHotspotDetectedInDecoy(contactMap, maxCap, hotspot3D.getPatients().size());
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
            hotspot3D.mergeHotspot(hotspots.iterator().next()); // add mutations
            hotspots3D.add(hotspot3D);
        }
        
        return Collections.singletonMap(protein, hotspots3D);
    }
    
    protected int[] getMutationCountsOnProtein(Map<Integer, Hotspot> mapResidueHotspot, int len) {
        int[] ret = new int[len+1];
        for (Map.Entry<Integer, Hotspot> entry : mapResidueHotspot.entrySet()) {
            ret[entry.getKey()] = entry.getValue().getPatients().size();
        }
        return ret;//Arrays.asList(ArrayUtils.toObject(ret));
    }
    
    protected int[][] generateDecoys(int[] counts, int left, int right, int times) {
        int[][] decoys = new int[times][];
        for (int i=0; i<times; i++) {
            decoys[i] = Arrays.copyOf(counts, counts.length);
            shuffleArray(decoys[i], left, right);
        }
        return decoys;
    }
    
    public static void shuffleArray(int[] a, int left, int right) {
        Random random = new Random();
        random.nextInt();
        for (int i = left; i < right; i++) {
          int change = i + random.nextInt(right - i);
          swap(a, i, change);
        }
    }

    private static void swap(int[] a, int i, int change) {
        int helper = a[i];
        a[i] = a[change];
        a[change] = helper;
    }
    
    static interface DetectedInDecoy {
        boolean isDetectedInDecoy(final int[] decoy);
    }
    
    private static class StructureHotspotDetectedInDecoy implements DetectedInDecoy {
        private final ContactMap contactMap;
        private final int maxCap;
        private final int targetCount;
        StructureHotspotDetectedInDecoy(final ContactMap contactMap, final int maxCap, final int targetCount) {
            this.contactMap = contactMap;
            this.maxCap = maxCap;
            this.targetCount = targetCount;
        }
        public boolean isDetectedInDecoy(final int[] decoy) {
            boolean[][] graph = contactMap.getContact();
            int l = contactMap.getProteinLeft();
            int r = contactMap.getProteinRight();
            for (int i=l; i<r; i++) {
                int count = 0;
                for (int j=l; j<r; j++) {
                    if (graph[i][j]) {
                        int c = decoy[j];
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
    }
    
    protected double getP(final DetectedInDecoy detectedInDecoy, int[][] decoyCountsList) {
        final AtomicInteger d = new AtomicInteger(0);
        
        int nDecoy = decoyCountsList.length;
        int nThread = 50;
        int bin = nDecoy / nThread;
        Thread[] threads = new Thread[nThread];
        for (int i = 0; i < nThread; i++) {
            final int[][] decoys = java.util.Arrays.copyOfRange(decoyCountsList, i*bin, (i+1)*bin);
            threads[i] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int[] decoy : decoys) {
                        if (detectedInDecoy.isDetectedInDecoy(decoy)) {
                            d.incrementAndGet();
                        }
                    }
                }
            });
            threads[i].start();
        }
        
        for (Thread thread : threads) {
            try {
                thread.join();
            }catch (InterruptedException ex) {
                throw new RuntimeException(ex);
            }
        }
        
        
        return 1.0 * d.get() / decoyCountsList.length;
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
    
    private Map<MutatedProtein3D,ContactMap> getContactMaps(MutatedProtein protein) throws HotspotException {
        try {
            Map<MutatedProtein3D,ContactMap> contactMaps = new HashMap<MutatedProtein3D,ContactMap>();
            Map<MutatedProtein3D,List<PdbUniprotAlignment>> mapAlignments = getPdbUniprotAlignments(protein);
            for (Map.Entry<MutatedProtein3D,List<PdbUniprotAlignment>> entryMapAlignments : mapAlignments.entrySet()) {
                MutatedProtein3D protein3D = entryMapAlignments.getKey();
                List<PdbUniprotAlignment> alignments = entryMapAlignments.getValue();
                OneToOneMap<Integer, Integer> pdbUniprotResidueMapping = getPdbUniprotResidueMapping(alignments);
                Map<Integer, Set<Integer>> pdbContactMap = getPdbContactMap(protein3D);
                
                ContactMap contactMap = getUniProtContactMap(pdbContactMap, pdbUniprotResidueMapping, protein.getProteinLength());
                
                contactMaps.put(protein3D, contactMap);
            }
            
            return contactMaps;
        } catch (DaoException e) {
            throw new HotspotException(e);
        }
    }
    
    private ContactMap getUniProtContactMap(Map<Integer, Set<Integer>> pdbContactMap, OneToOneMap<Integer, Integer> pdbUniprotResidueMapping, int proteinLength) {
        ContactMap contactMap = new ContactMap(proteinLength);
        contactMap.setProteinLeft(pdbUniprotResidueMapping.getSmallestValue());
        contactMap.setProteinRight(pdbUniprotResidueMapping.getLargestValue());
        
        boolean[][] matrix = contactMap.getContact();
        
        for (Map.Entry<Integer, Set<Integer>> entryPdbContactMap : pdbContactMap.entrySet()) {
            int pdbPos = entryPdbContactMap.getKey();
            if (!pdbUniprotResidueMapping.hasKey(pdbPos)) {
                continue;
            }
            
            int uniprotPos = pdbUniprotResidueMapping.getByKey(pdbPos);
            if (uniprotPos > proteinLength) {
                System.err.println("Mapped Protein resisue longer than protein length");
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
                if (uniprotNeighbor > proteinLength) {
                    System.err.println("Mapped Protein resisue longer than protein length");
                    continue;
                }
                matrix[uniprotPos][uniprotNeighbor] = true;
                matrix[uniprotNeighbor][uniprotPos] = true;
            }
        }
        
        return contactMap;
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
    
    final class OneToOneMap<K extends Comparable, V extends Comparable> {
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
        
        K getSmallestKey() {
            return (K) Collections.min(keyToValMap.keySet());
        }
        
        K getLargestKey() {
            return (K) Collections.max(keyToValMap.keySet());
        }
        
        V getSmallestValue() {
            return (V) Collections.min(valToKeyMap.keySet());
        }
        
        V getLargestValue() {
            return (V) Collections.max(valToKeyMap.keySet());
        }
    }
    
    final class ContactMap {
        private boolean[][] contact;
        private int proteinLeft, proteinRight;
        
        ContactMap(int len) {
            contact = new boolean[len+1][len+1];
        }

        public boolean[][] getContact() {
            return contact;
        }

        public int getProteinLeft() {
            return proteinLeft;
        }

        public void setProteinLeft(int proteinLeft) {
            this.proteinLeft = proteinLeft;
        }

        public int getProteinRight() {
            return proteinRight;
        }

        public void setProteinRight(int proteinRight) {
            this.proteinRight = proteinRight;
        }

        
    }
}
