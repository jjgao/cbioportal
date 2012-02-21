package org.mskcc.portal.mut_diagram.servlet;

import com.google.common.collect.ImmutableList;
import org.apache.log4j.Logger;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.map.type.CollectionType;
import org.codehaus.jackson.map.type.TypeFactory;
import org.mskcc.cgds.dao.DaoException;
import org.mskcc.cgds.dao.DaoGeneOptimized;
import org.mskcc.cgds.model.ExtendedMutation;
import org.mskcc.portal.mut_diagram.*;
import org.mskcc.portal.mut_diagram.impl.CacheFeatureService;
import org.mskcc.portal.mut_diagram.impl.CgdsIdMappingService;
import org.mskcc.portal.mut_diagram.impl.PfamGraphicsCacheLoader;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;

import static com.google.common.collect.Lists.newArrayList;
import static com.google.common.collect.Maps.newHashMap;
import static org.codehaus.jackson.map.DeserializationConfig.Feature.ACCEPT_SINGLE_VALUE_AS_ARRAY;

/**
 * Mutation diagram data servlet.
 */
public final class OncoVizDataServlet extends HttpServlet {
    private static final Logger logger = Logger.getLogger(OncoVizDataServlet.class);
    /** Default serial version UID. */
    private static final long serialVersionUID = 1L;
    private static final List<Sequence> EMPTY = Collections.emptyList();

    private final ObjectMapper objectMapper;
    private final FeatureService featureService;
    private final IdMappingService idMappingService;

    public OncoVizDataServlet() {
        objectMapper = new ObjectMapper();
        objectMapper.configure(ACCEPT_SINGLE_VALUE_AS_ARRAY, true);

        PfamGraphicsCacheLoader cacheLoader = new PfamGraphicsCacheLoader(objectMapper);
        featureService = new CacheFeatureService(cacheLoader);

        try {
            idMappingService = new CgdsIdMappingService(DaoGeneOptimized.getInstance());
        }
        catch (DaoException e) {
            throw new RuntimeException("could not create id mapping service", e);
        }
    }

    @Override
    protected void doPost(final HttpServletRequest request, final HttpServletResponse response) 
            throws ServletException, IOException {
        String hugoGeneSymbol = request.getParameter("hugoGeneSymbol");

        List<String> uniProtIds = idMappingService.getUniProtIds(hugoGeneSymbol);
        if (uniProtIds.isEmpty()) {
            writeSequencesToResponse(hugoGeneSymbol, EMPTY, response);
            return;
        }

        String uniProtId = uniProtIds.get(0);
        List<Sequence> sequences = featureService.getFeatures(uniProtId);
        if (sequences.isEmpty()) {
            writeSequencesToResponse(hugoGeneSymbol, EMPTY, response);
            return;
        }

        Sequence sequence = sequences.get(0);
        if (sequence.getMetadata() == null) {
            Map<String, Object> metadata = newHashMap();
            sequence.setMetadata(metadata);
        }
        sequence.getMetadata().put("hugoGeneSymbol", hugoGeneSymbol);
        sequence.getMetadata().put("uniProtId", uniProtId);
 
        List<ExtendedMutation> mutations = readMutations(request.getParameter("mutations"));

        List<Markup> markups = createMarkups(mutations);
        writeSequencesToResponse(hugoGeneSymbol, ImmutableList.of(sequence.withMarkups(markups)), response);
    }

    /**
     * Creats Markups From Specified Mutations.
     */
    private List<Markup> createMarkups(List<ExtendedMutation> mutations) {
        List<Markup> markups = newArrayList();
        HashMap<Integer, Integer> locationMap = new HashMap<Integer, Integer>();
        for (ExtendedMutation mutation: mutations) {
            String aaChange = mutation.getAminoAcidChange();
            try {
                int location = Integer.valueOf(aaChange.replaceAll("[A-Za-z\\.*]+", ""));
                int yoffset = 1;
                if (locationMap.containsKey(location)) {
                    yoffset = locationMap.get(location) + 1;
                }
                locationMap.put(location, yoffset);
                Markup markup = new Markup();
                markup.setDisplay("true");
                markup.setStart(location);
                markup.setEnd(location);
                markup.setColour(ImmutableList.of(mutation.getColor()));
                markup.setLineColour("#babdb6");
                markup.setHeadStyle("diamond");
                markup.setYoffset(yoffset);

                markup.setV_align("top");
                markup.setType("mutation");
                markup.setMetadata(new HashMap<String, Object>());
                markup.getMetadata().put("label", mutation.getAminoAcidChange());
                markups.add(markup);
            }
            catch (NumberFormatException e) {
                logger.warn("ignoring extended mutation " + aaChange + ", no location information");
            }
        }
        return markups;
    }

    private void writeSequencesToResponse(String geneSymbol,
              final List<Sequence> sequences, final HttpServletResponse response)
            throws IOException {
        response.setContentType("text/html");
        PrintWriter writer = response.getWriter();
        writer.write("<html>");
        writer.write("<body>");
        writer.write("<h2>OncoViz</h2>");
        writer.write("Almost there...<P>");
        writer.write("<form action='onco_viz.do' method='POST'>");
        writer.write("<input type='hidden' name='hugoGeneSymbol' value='"
            + geneSymbol + "'>");
        writer.write("<textarea name=\"json\" rows=\"20\" cols=\"60\">");
        StringWriter strWriter = new StringWriter();
        objectMapper.writeValue(strWriter, sequences);
        writer.write(strWriter.toString());
        writer.write("</textarea>");
        writer.write("<P><input type=\"submit\">");
        writer.write("</form>");
        writer.write("</body>");
        writer.write("</html>");
    }

    /**
     * Reads in mutations from the Request Object.
     */
    List<ExtendedMutation> readMutations(final String value) {
        List<ExtendedMutation> mutations = new ArrayList<ExtendedMutation>();
        if (value != null) {
            String lines[] = value.split("\n");
            for (String line: lines) {
                if (line.trim().length()>0) {
                    //  Split by any white space
                    String parts[] = line.split("\\s+");

                    String caseId = parts[0];
                    String aaChange = parts[1];
                    String color = "BLACK";
                    try {
                        color = parts[2];
                    } catch (ArrayIndexOutOfBoundsException e) {
                    }
                    ExtendedMutation mutation = new ExtendedMutation();
                    mutation.setAminoAcidChange(aaChange);
                    mutation.setCaseId(caseId);
                    mutation.setColor(color);
                    mutations.add(mutation);
                }
            }
        }
        return mutations;
    }
}
