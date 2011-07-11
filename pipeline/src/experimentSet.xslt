<?xml version="1.0" encoding="UTF-8"?>

<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY amp "&#64257;">
]>




<xsl:stylesheet 

	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 

	xmlns:xs="http://www.w3.org/2001/XMLSchema"
	xmlns:xalan="http://xml.apache.org/xalan"
	xmlns:redirect="http://xml.apache.org/xalan/redirect" 
	extension-element-prefixes="redirect" 
	exclude-result-prefixes="xs xalan"
	version="1.0">
	
	<xsl:output method="text"/>
	<xsl:param name="checkDir" select="''"/>
	<xsl:param name="generate_params_feature" select="'false'"/>
	<xsl:param name="generate_params_behavior" select="'false'"/>
	<xsl:param name="generate_scripts" select="'false'"/>
	
	<xsl:variable name="schema" select="document('../../../data/experimentSet.xsd')" />
	<xsl:variable name="document" select="."/>
	<xsl:variable name="outputDir" select="$document/experimentSet/@dir"/>
	<xsl:variable name="diffFile" select="concat($checkDir, '/', 'diff.sh')"/>

<!--	
	<xsl:variable name="BehaviorFilename" select="concat(Params/@dir,'/',Params/BehaviorParams/@file)"/>
	<xsl:variable name="FeatureFilename" select="concat(Params/@dir,'/',Params/FeatureParams/@file)"/>
-->	
	
	
	<xsl:template name="valueOrDefault">
		<xsl:param name="value"/>
		<xsl:param name="default"/>
		
		<xsl:choose>
			<xsl:when test="$value">
				<xsl:value-of select="$value" />
			</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="$default" />
			</xsl:otherwise>
		</xsl:choose>		
	</xsl:template>
	
	 <xsl:template match="experimentSet">	 
		<xsl:choose>
			<xsl:when test="$checkDir=''">
				<xsl:call-template name="processExperimentSet">
					<xsl:with-param name="outputDir" select="@dir"/>
  				</xsl:call-template>
			</xsl:when>
			<xsl:otherwise>
				<redirect:write select="$diffFile">
					<xsl:text># all diffs between checkDir and original outputDir, this | wc -l should be 0
</xsl:text>
				</redirect:write>
				<xsl:call-template name="processExperimentSet">
					<xsl:with-param name="outputDir" select="$checkDir"/>
  				</xsl:call-template>
			
			</xsl:otherwise>
		</xsl:choose>
	</xsl:template>
	
	<xsl:template name="processExperimentSet">
		<xsl:param name="outputDir"/>
		
		<xsl:call-template name="generateFile">
			<xsl:with-param name="filename">
				<xsl:call-template name="getBehaviorFilename"><xsl:with-param name="outputDir" select="$outputDir"/></xsl:call-template>
			</xsl:with-param>
			<xsl:with-param name="generate" select="$generate_params_behavior"/>
			<xsl:with-param name="func" select="'generateBehaviorParams'"/>
			<xsl:with-param name="display" select="1"/>
		</xsl:call-template>
		
		<xsl:call-template name="generateFile">
			<xsl:with-param name="filename">
				<xsl:call-template name="getFeatureFilename"><xsl:with-param name="outputDir" select="$outputDir"/></xsl:call-template>
			</xsl:with-param>			
			<xsl:with-param name="generate" select="$generate_params_feature"/>
			<xsl:with-param name="func" select="'generateFeatureParams'"/>
			<xsl:with-param name="display" select="1"/>
		</xsl:call-template>
		
		<xsl:call-template name="generateFile">
			<xsl:with-param name="filename" select="concat($outputDir,'/','master.sh')"/>
			<xsl:with-param name="generate" select="$generate_scripts"/>
			<xsl:with-param name="func" select="'generateMasterScript'"/>
			<xsl:with-param name="display" select="1"/>
			<xsl:with-param name="funcparam" select="$outputDir"/>
		</xsl:call-template>

		
		<xsl:for-each select="//Model">
			<xsl:call-template name="generateFile">
				<xsl:with-param name="filename" select="concat($outputDir,'/',@dir,'/','process.sh')"/>
				<xsl:with-param name="generate" select="$generate_scripts"/>
				<xsl:with-param name="func" select="'generateExperimentScripts'"/>
				<xsl:with-param name="display" select="1"/>				
				<xsl:with-param name="funcparam" select="$outputDir"/>
			</xsl:call-template>
		</xsl:for-each>
		
	</xsl:template>

	<xsl:template name="appendDiff">
	<xsl:param name="filename"/>
		<xsl:variable name="origName" select="concat($outputDir, substring-after($filename, $checkDir))"/>
		<redirect:write append="yes" select="$diffFile">
			<xsl:text>diff </xsl:text>
			<xsl:value-of select="$filename"/>
			<xsl:text> </xsl:text>
			<xsl:value-of select="$origName"/>
			<xsl:text>
</xsl:text>			
		</redirect:write>
	</xsl:template>
	
	
	<xsl:template name="generateFile">
	<xsl:param name="filename"/>
	<xsl:param name="generate"/>
	<xsl:param name="display" select="1"/>
	<xsl:param name="func"/>
	<xsl:param name="funcparam" select="0"/>	
						<xsl:call-template name="appendDiff"><xsl:with-param name="filename" select="$filename"></xsl:with-param></xsl:call-template>

		<xsl:choose>
			<xsl:when test="$generate = 'true'">
			
			
				<xsl:if test="$checkDir != ''">
					<xsl:call-template name="appendDiff"><xsl:with-param name="filename" select="$filename"></xsl:with-param></xsl:call-template>
				</xsl:if>
			
				<xsl:text>
writing file </xsl:text><xsl:value-of select="$filename"/>
			              <redirect:write select="$behaviorFilename">
					<xsl:call-template name="generateFileContent">
						<xsl:with-param name="func" select="$func"/>
						<xsl:with-param name="funcparam" select="$funcparam"/>
					</xsl:call-template>
		   	             </redirect:write>			   	             
				<xsl:text>
				done.</xsl:text>
			</xsl:when>
			<xsl:otherwise>
				<xsl:if test="$display">
<xsl:text>
### Content </xsl:text><xsl:value-of select="$filename"/><xsl:text> BEGIN ###

</xsl:text>
					<xsl:call-template name="generateFileContent">
						<xsl:with-param name="func" select="$func"/>
						<xsl:with-param name="funcparam" select="$funcparam"/>
					</xsl:call-template>
<xsl:text>
### Content </xsl:text><xsl:value-of select="$filename"/><xsl:text> END ###

</xsl:text>
				</xsl:if>
			</xsl:otherwise>
		</xsl:choose>	
	</xsl:template>
	
	<xsl:template name="generateFileContent">
	<xsl:param name="func"/>
	<xsl:param name="funcparam"/>	
		<xsl:choose>
			<xsl:when test="$func='generateFeatureParams'">
				<xsl:call-template name="generateFeatureParams"/>
			</xsl:when>
			<xsl:when test="$func='generateBehaviorParams'">
				<xsl:call-template name="generateBehaviorParams"/>
			</xsl:when>
			<xsl:when test="$func='generateMasterScript'">
				<xsl:call-template name="generateMasterScript">
					<xsl:with-param name="outputDir" select="$funcparam"/>
				</xsl:call-template>
			</xsl:when>
			<xsl:when test="$func='generateExperimentScripts'">
				<xsl:call-template name="generateExperimentScripts">
					<xsl:with-param name="outputDir" select="$funcparam"/>
				</xsl:call-template>				
			</xsl:when>
			<xsl:when test="$func='generateFileList'">
				<xsl:call-template name="generateFileList">
					<xsl:with-param name="FileSet" select="$funcparam"/>
				</xsl:call-template>
			</xsl:when>
			<xsl:when test="$func='writeContent'">
				<xsl:call-template name="writeContent">
					<xsl:with-param name="content" select="$funcparam"/>
				</xsl:call-template>
			</xsl:when>
			<xsl:otherwise>
				Error: This should never happen: No func specified!
			</xsl:otherwise>
		</xsl:choose>
	</xsl:template>
	
		
	<xsl:template name="generateFeatureParams">
		<xsl:for-each select="Params/FeatureParams//Feature">
			<xsl:variable name="currFeature" select="."/>
			
			<xsl:value-of select="$currFeature/@name"/><xsl:text>: </xsl:text>

			<xsl:for-each select="$schema/xs:schema/xs:attributeGroup[@name='boutLevelFeatures']//xs:attribute">

				<xsl:value-of select="@name"/>
				<xsl:variable name="currBoutFeatureName" select="@name"/>
				<xsl:text>=</xsl:text>
				<xsl:call-template  name="valueOrDefault"><xsl:with-param name="value" select="$currFeature/@*[local-name(.)=$currBoutFeatureName]"/><xsl:with-param name="default" select="@default"/></xsl:call-template>				
				
				<xsl:if test="position() &lt; last()">, </xsl:if>
			</xsl:for-each>

			<xsl:text>
</xsl:text>
		</xsl:for-each>
	</xsl:template>
	

	<xsl:template name="generateBehaviorParams">
		<xsl:text>9	1
</xsl:text>
		<xsl:for-each select="Params/BehaviorParams//BehaviorGroup">
			<xsl:text>*</xsl:text>
			<xsl:value-of select="@name"/>
			<xsl:text>	</xsl:text>
			<xsl:call-template  name="valueOrDefault"><xsl:with-param name="value" select="@classifierMethod"/><xsl:with-param name="default" select="$schema/xs:schema/xs:complexType[@name='BehaviorGroupType']/xs:attribute[@name='classifierMethod']/@default"/></xsl:call-template>
			<xsl:text>	</xsl:text>
			<xsl:call-template  name="valueOrDefault"><xsl:with-param name="value" select="@isMultiClass"/><xsl:with-param name="default" select="$schema/xs:schema/xs:complexType[@name='BehaviorGroupType']/xs:attribute[@name='isMultiClass']/@default"/></xsl:call-template>
			<xsl:text>
</xsl:text>

			<xsl:for-each select="//BehaviorClass">
				<xsl:value-of select="@name"/>
				<xsl:text>	</xsl:text>
				<xsl:value-of select="@name"/>
				<xsl:text>	</xsl:text>
				<xsl:call-template  name="valueOrDefault"><xsl:with-param name="value" select="@classifierMethod"/><xsl:with-param name="default" select="$schema/xs:schema/xs:complexType[@name='BehaviorClassType']/xs:attribute[@name='classifierMethod']/@default"/></xsl:call-template>
				<xsl:text>	</xsl:text>
				<xsl:value-of select="@color"/>
				<xsl:text>
</xsl:text>
			</xsl:for-each>
		</xsl:for-each>
	</xsl:template>

	<xsl:template name="generateQsubCmd">
	<xsl:param name="scriptName"/>
			<xsl:text>qsub QSUB PARAMS </xsl:text><xsl:value-of  select="$scriptName"/><xsl:text>
</xsl:text>
	</xsl:template>

	
	<xsl:template name="generateMasterScript">
	<xsl:param name="outputDir" select="'**UNDEFINED**'"/>
		
		<xsl:for-each select="//Model">
			<xsl:call-template name="generateQsubCmd"><xsl:with-param name="scriptName" select="concat($outputDir,'/',@dir,'/','process.sh')"/></xsl:call-template>			
		</xsl:for-each>
	</xsl:template>
	
	
	<xsl:template name="getBehaviorFilename">
	<xsl:param name="outputDir" select="'****UNDEFINED****'"/>
		<xsl:value-of select="concat($outputDir,'/',$document/experimentSet//Params/@dir,'/',$document/experimentSet/Params/BehaviorParams/@file)"/>	
	</xsl:template>
	<xsl:template name="getFeatureFilename">
	<xsl:param name="outputDir" select="'****UNDEFINED****'"/>
		<xsl:value-of select="concat($outputDir,'/',$document/experimentSet//Params/@dir,'/',$document/experimentSet/Params/FeatureParams/@file)"/>	
	</xsl:template>
	
	<xsl:template name="generateFileList">
	<xsl:param name="FileSet"/>
		<xsl:for-each select="$FileSet">
			<xsl:value-of select="@labelFile"/>
			<xsl:text>
</xsl:text>
		</xsl:for-each>
	</xsl:template>

	<xsl:template name="generateTrainInvocation">
	<xsl:param name="outputDir" select="'***UNDEFINED***'"/>
	<xsl:param name="ModelDir"/>
	<xsl:param name="debugDir"/>
	<xsl:param name="debugParams"/>
	<xsl:param name="ModelFile"/>
	<xsl:param name="TrainFile"/>
	<xsl:param name="TrainSet"/>

		<xsl:variable name="destTrainFile" select="concat($outputDir, '/', $ModelDir, '/', $TrainFile)"/>

		<xsl:if test="$TrainSet">
			<xsl:call-template name="generateFile">
				<xsl:with-param name="filename" select="$destTrainFile"/>
				<xsl:with-param name="generate" select="$generate_scripts"/>
				<xsl:with-param name="func" select="'generateFileList'"/>
				<xsl:with-param name="display" select="1"/>
				<xsl:with-param name="funcparam" select="$TrainSet"/>
			</xsl:call-template>
		</xsl:if>	

		<xsl:text>./svm_struct_train</xsl:text>
		<xsl:text> -F </xsl:text><xsl:call-template name="getFeatureFilename"><xsl:with-param name="outputDir" select="$outputDir"/></xsl:call-template>
		<xsl:text> -B </xsl:text><xsl:call-template name="getBehaviorFilename"><xsl:with-param name="outputDir" select="$outputDir"/></xsl:call-template>
		<xsl:text> -D</xsl:text><xsl:value-of select="$debugParams"/><xsl:text> </xsl:text><xsl:value-of select="$debugDir"/>
		<xsl:text> </xsl:text><xsl:value-of select="$destTrainFile"/>
		<xsl:text> </xsl:text><xsl:value-of select="concat($outputDir, '/', $ModelDir, '/', $ModelFile)"/>
		<xsl:text>
</xsl:text>
	</xsl:template>

	
	<xsl:template name="generateTestInvocation">
	<xsl:param name="outputDir" select="'***UNDEFINED***'"/>
	<xsl:param name="ModelDir"/>
	<xsl:param name="debugDir"/>
	<xsl:param name="debugParams"/>
	<xsl:param name="ModelFile"/>
	<xsl:param name="file"/>
	<xsl:param name="FileSet"/>

		<xsl:variable name="destFile" select="concat($outputDir, '/', $ModelDir, '/', $file)"/>

		<xsl:if test="$FileSet">
			<xsl:call-template name="generateFile">
				<xsl:with-param name="filename" select="$destFile"/>
				<xsl:with-param name="generate" select="$generate_scripts"/>
				<xsl:with-param name="func" select="'generateFileList'"/>
				<xsl:with-param name="display" select="1"/>
				<xsl:with-param name="funcparam" select="$FileSet"/>
			</xsl:call-template>
		</xsl:if>	

		<xsl:text>./svm_struct_classify -c 1000 -e 5</xsl:text>
		<xsl:text> -F </xsl:text><xsl:call-template name="getFeatureFilename"><xsl:with-param name="outputDir" select="$outputDir"/></xsl:call-template>
		<xsl:text> -B </xsl:text><xsl:call-template name="getBehaviorFilename"><xsl:with-param name="outputDir" select="$outputDir"/></xsl:call-template>
<!--		<xsl:text> -D</xsl:text><xsl:value-of select="$debugParams"/><xsl:text> </xsl:text><xsl:value-of select="$debugDir"/>-->
		<xsl:text> -O </xsl:text><xsl:value-of select="$destFile"/>
		<xsl:text> </xsl:text><xsl:value-of select="concat($outputDir, '/', $ModelDir, '/', $ModelFile)"/>
		<xsl:text>
</xsl:text>
	</xsl:template>
	
	<xsl:template name="writeContent">
	<xsl:param name="content"/>
		<xsl:value-of select="$content"/>
	</xsl:template>	
	
	<xsl:template name="generateExperimentScripts">
	<xsl:param name="outputDir" select="'**UNDEFINED**'"/>
	
		<xsl:variable name="ModelDir" select="@dir"/>
		<xsl:choose>
			<xsl:when test="count(TrainList) > 0">
				<xsl:call-template name="generateTrainInvocation">
					<xsl:with-param name="outputDir" select="$outputDir"/>
					<xsl:with-param name="ModelDir" select="$ModelDir"/>
					<xsl:with-param name="debugDir" select="concat($outputDir, '/', $ModelDir, '/', TrainList/@debugDir)"/>
					<xsl:with-param name="debugParams" select="TrainList/@debugParams"/>
					<xsl:with-param name="ModelFile" select="'LearnedModel.txt'"/>
					<xsl:with-param name="TrainFile" select="'Trainlist.txt'"/>
					<xsl:with-param name="TrainSet" select="TrainList/TrainItem"/>
				</xsl:call-template>
				
				<xsl:if test="TrainList/@resultLabelDir">
					<xsl:call-template name="generateTestInvocation">
						<xsl:with-param name="outputDir" select="$outputDir"/>
						<xsl:with-param name="ModelDir" select="$ModelDir"/>
						<xsl:with-param name="debugDir" select="concat($outputDir, '/', $ModelDir, '/', TrainList/@debugDir)"/>
						<xsl:with-param name="debugParams" select="TrainList/@debugParams"/>
						<xsl:with-param name="ModelFile" select="'LearnedModel.txt'"/>
						<xsl:with-param name="file" select="'Trainlist.txt'"/>
					</xsl:call-template>
				</xsl:if>												

			</xsl:when>
			<xsl:when test="count(CrossValidationList) > 0">
				<xsl:variable name="CrossValList" select="CrossValidationList"/>
				
				<xsl:for-each select="CrossValidationList/Item">
					<xsl:variable name="currPos" select="position()"/>

					<xsl:variable name="processThis" select="concat($outputDir, '/', $ModelDir, '/', 'process', '.', $currPos, '.sh')"/>
					<xsl:call-template name="generateQsubCmd"><xsl:with-param name="scriptName" select="$processThis"/></xsl:call-template>

					
<xsl:call-template name="generateFile">
						<xsl:with-param name="filename" select="$processThis"/>
						<xsl:with-param name="generate" select="$generate_scripts"/>
						<xsl:with-param name="func" select="'writeContent'"/>
						<xsl:with-param name="display" select="1"/>
						
<xsl:with-param name="funcparam">
		
					<xsl:call-template name="generateTrainInvocation">
								<xsl:with-param name="outputDir" select="$outputDir"/>
								<xsl:with-param name="ModelDir" select="$ModelDir"/>
								<xsl:with-param name="debugDir" select="concat($outputDir, '/', $ModelDir, '/', $CrossValList/@debugDir)"/>
								<xsl:with-param name="debugParams" select="$CrossValList/@debugParams"/>
								<xsl:with-param name="ModelFile" select="concat('LearnedModel.txt','.',$currPos)"/>
								<xsl:with-param name="TrainFile" select="concat('Trainlist.txt', '.', $currPos)"/>
								<xsl:with-param name="TrainSet" select="$CrossValList/Item[position() != $currPos]"/>
							</xsl:call-template>

							<xsl:call-template name="generateTestInvocation">
								<xsl:with-param name="outputDir" select="$outputDir"/>
								<xsl:with-param name="ModelDir" select="$ModelDir"/>
								<xsl:with-param name="debugDir" select="concat($outputDir, '/', $ModelDir, '/', $CrossValList/@debugDir)"/>
								<xsl:with-param name="debugParams" select="$CrossValList/@debugParams"/>
								<xsl:with-param name="ModelFile" select="concat('LearnedModel.txt','.',$currPos)"/>
								<xsl:with-param name="file" select="concat('Trainlist.txt', '.', $currPos)"/>
								<xsl:with-param name="FileSet" select="$CrossValList/Item[position() = $currPos]"/>
							</xsl:call-template>
					
						
</xsl:with-param>
					
</xsl:call-template>
				</xsl:for-each>
				
			</xsl:when>
			<xsl:otherwise>
				Error: This should never happen: Input file violates schema!
			</xsl:otherwise>
		</xsl:choose>
		<xsl:for-each select="TestList">
		
					<xsl:call-template name="generateTestInvocation">
						<xsl:with-param name="outputDir" select="$outputDir"/>
						<xsl:with-param name="ModelDir" select="$ModelDir"/>
<!--						<xsl:with-param name="debugDir" select="concat($outputDir, '/', @dir, '/', TrainList/@debugDir)"/>->
<!-						<xsl:with-param name="debugParams" select="TrainList/@debugParams"/>-->
						<xsl:with-param name="ModelFile" select="'LearnedModel.txt'"/>
						<xsl:with-param name="file" select="concat(@resultLabelDir, '/', 'TestList.txt')"/>
						<xsl:with-param name="FileSet" select="TestItem"/>
					</xsl:call-template>

		</xsl:for-each>
	</xsl:template>

</xsl:stylesheet>
