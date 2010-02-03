<?xml version='1.0' encoding='UTF-8'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version='2.0'
                xmlns="http://www.w3.org/1999/xhtml">
<xsl:output 
      doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN"
      doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"
      method='xml' encoding='UTF-8' indent='yes'/>

<xsl:template match='QPCR'>
  <html lang="en">

    <head>
      <meta http-equiv="Content-Type" content="text/html;charset=UTF-8" />
      <title> qPCR results </title>
      <link rel="stylesheet" href="pyqpcr.css" type="text/css"/>
    </head>

    <body>
      <h1> qPCR results </h1>
      <p><br/><br/></p>
      <xsl:apply-templates select='PLATE'/>
    </body>

  </html>
</xsl:template>

<xsl:template match='PLATE'>
  <h2 > <xsl:value-of select="@NAME"/> </h2>
  <table>
      <tr>
        <th class="th_title">Well</th>
        <th class="th_title">Type</th>
        <th class="th_title">Target</th>
        <th class="th_title">Sample</th>
        <th class="th_title">Ct</th>
        <th class="th_title">Ctmean</th>
        <th class="th_title">Ctdev</th>
        <th class="th_title">Amount</th>
        <th class="th_title">Efficiency</th>
        <th class="th_title">NRQ</th>
        <th class="th_title">NRQerror</th>
      </tr>
      <xsl:apply-templates select='WELL'/>
  </table>
  <p><br/></p>
</xsl:template>
 
<xsl:template match='WELL'>
  <tr>
    <th class="th_inside"> <xsl:value-of select="NAME"/></th>
    <th class="th_inside"> <xsl:value-of select="TYPE"/></th>
    <th class="th_inside"> <xsl:value-of select="TARGET"/></th>
    <th class="th_inside"> <xsl:value-of select="SAMPLE"/></th>
    <th class="th_inside"> <xsl:value-of select="@CT"/></th>
    <th class="th_inside"> <xsl:value-of select="@CTMEAN"/></th>
    <th class="th_inside"> <xsl:value-of select="@CTDEV"/></th>
    <th class="th_inside"> <xsl:value-of select="@AMOUNT"/></th>
    <th class="th_inside"> <xsl:apply-templates select='TARGET'/></th>
    <th class="th_inside"> <xsl:value-of select="@NRQ"/></th>
    <th class="th_inside"> <xsl:value-of select="@NRQERROR"/></th>
  </tr>
</xsl:template>

<xsl:template match='TARGET'>
  <xsl:value-of select="@EFF"/>&#177;<xsl:value-of select="@PM"/>
</xsl:template>

</xsl:stylesheet>
