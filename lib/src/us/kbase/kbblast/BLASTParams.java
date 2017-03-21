
package us.kbase.kbblast;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: BLAST_Params</p>
 * <pre>
 * BLAST Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_one_sequence",
    "input_one_ref",
    "input_many_ref",
    "input_msa_ref",
    "output_one_name",
    "output_filtered_name",
    "ident_thresh",
    "e_value",
    "bitscore",
    "overlap_fraction",
    "maxaccepts",
    "output_extra_format",
    "rounds"
})
public class BLASTParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_one_sequence")
    private String inputOneSequence;
    @JsonProperty("input_one_ref")
    private String inputOneRef;
    @JsonProperty("input_many_ref")
    private String inputManyRef;
    @JsonProperty("input_msa_ref")
    private String inputMsaRef;
    @JsonProperty("output_one_name")
    private String outputOneName;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("ident_thresh")
    private Double identThresh;
    @JsonProperty("e_value")
    private Double eValue;
    @JsonProperty("bitscore")
    private Double bitscore;
    @JsonProperty("overlap_fraction")
    private Double overlapFraction;
    @JsonProperty("maxaccepts")
    private Double maxaccepts;
    @JsonProperty("output_extra_format")
    private String outputExtraFormat;
    @JsonProperty("rounds")
    private Double rounds;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public BLASTParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_one_sequence")
    public String getInputOneSequence() {
        return inputOneSequence;
    }

    @JsonProperty("input_one_sequence")
    public void setInputOneSequence(String inputOneSequence) {
        this.inputOneSequence = inputOneSequence;
    }

    public BLASTParams withInputOneSequence(String inputOneSequence) {
        this.inputOneSequence = inputOneSequence;
        return this;
    }

    @JsonProperty("input_one_ref")
    public String getInputOneRef() {
        return inputOneRef;
    }

    @JsonProperty("input_one_ref")
    public void setInputOneRef(String inputOneRef) {
        this.inputOneRef = inputOneRef;
    }

    public BLASTParams withInputOneRef(String inputOneRef) {
        this.inputOneRef = inputOneRef;
        return this;
    }

    @JsonProperty("input_many_ref")
    public String getInputManyRef() {
        return inputManyRef;
    }

    @JsonProperty("input_many_ref")
    public void setInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
    }

    public BLASTParams withInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
        return this;
    }

    @JsonProperty("input_msa_ref")
    public String getInputMsaRef() {
        return inputMsaRef;
    }

    @JsonProperty("input_msa_ref")
    public void setInputMsaRef(String inputMsaRef) {
        this.inputMsaRef = inputMsaRef;
    }

    public BLASTParams withInputMsaRef(String inputMsaRef) {
        this.inputMsaRef = inputMsaRef;
        return this;
    }

    @JsonProperty("output_one_name")
    public String getOutputOneName() {
        return outputOneName;
    }

    @JsonProperty("output_one_name")
    public void setOutputOneName(String outputOneName) {
        this.outputOneName = outputOneName;
    }

    public BLASTParams withOutputOneName(String outputOneName) {
        this.outputOneName = outputOneName;
        return this;
    }

    @JsonProperty("output_filtered_name")
    public String getOutputFilteredName() {
        return outputFilteredName;
    }

    @JsonProperty("output_filtered_name")
    public void setOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
    }

    public BLASTParams withOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
        return this;
    }

    @JsonProperty("ident_thresh")
    public Double getIdentThresh() {
        return identThresh;
    }

    @JsonProperty("ident_thresh")
    public void setIdentThresh(Double identThresh) {
        this.identThresh = identThresh;
    }

    public BLASTParams withIdentThresh(Double identThresh) {
        this.identThresh = identThresh;
        return this;
    }

    @JsonProperty("e_value")
    public Double getEValue() {
        return eValue;
    }

    @JsonProperty("e_value")
    public void setEValue(Double eValue) {
        this.eValue = eValue;
    }

    public BLASTParams withEValue(Double eValue) {
        this.eValue = eValue;
        return this;
    }

    @JsonProperty("bitscore")
    public Double getBitscore() {
        return bitscore;
    }

    @JsonProperty("bitscore")
    public void setBitscore(Double bitscore) {
        this.bitscore = bitscore;
    }

    public BLASTParams withBitscore(Double bitscore) {
        this.bitscore = bitscore;
        return this;
    }

    @JsonProperty("overlap_fraction")
    public Double getOverlapFraction() {
        return overlapFraction;
    }

    @JsonProperty("overlap_fraction")
    public void setOverlapFraction(Double overlapFraction) {
        this.overlapFraction = overlapFraction;
    }

    public BLASTParams withOverlapFraction(Double overlapFraction) {
        this.overlapFraction = overlapFraction;
        return this;
    }

    @JsonProperty("maxaccepts")
    public Double getMaxaccepts() {
        return maxaccepts;
    }

    @JsonProperty("maxaccepts")
    public void setMaxaccepts(Double maxaccepts) {
        this.maxaccepts = maxaccepts;
    }

    public BLASTParams withMaxaccepts(Double maxaccepts) {
        this.maxaccepts = maxaccepts;
        return this;
    }

    @JsonProperty("output_extra_format")
    public String getOutputExtraFormat() {
        return outputExtraFormat;
    }

    @JsonProperty("output_extra_format")
    public void setOutputExtraFormat(String outputExtraFormat) {
        this.outputExtraFormat = outputExtraFormat;
    }

    public BLASTParams withOutputExtraFormat(String outputExtraFormat) {
        this.outputExtraFormat = outputExtraFormat;
        return this;
    }

    @JsonProperty("rounds")
    public Double getRounds() {
        return rounds;
    }

    @JsonProperty("rounds")
    public void setRounds(Double rounds) {
        this.rounds = rounds;
    }

    public BLASTParams withRounds(Double rounds) {
        this.rounds = rounds;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((((((((((((("BLASTParams"+" [workspaceName=")+ workspaceName)+", inputOneSequence=")+ inputOneSequence)+", inputOneRef=")+ inputOneRef)+", inputManyRef=")+ inputManyRef)+", inputMsaRef=")+ inputMsaRef)+", outputOneName=")+ outputOneName)+", outputFilteredName=")+ outputFilteredName)+", identThresh=")+ identThresh)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", overlapFraction=")+ overlapFraction)+", maxaccepts=")+ maxaccepts)+", outputExtraFormat=")+ outputExtraFormat)+", rounds=")+ rounds)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
