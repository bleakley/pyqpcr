<?xml version='1.0' encoding='UTF-8'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version='2.0'
                xmlns="http://www.w3.org/1999/xhtml">
<xsl:output 
      doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"
      doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"
      method='html' indent="yes" encoding='UTF-8'/>

<xsl:template match='QPCR'>
  <html lang="en">

    <head>
      <title> qPCR results </title>
    </head>

    <body>
      <table width="100%">
      <tr>
      <td align="center" bgcolor="#ff6600"><font size="5"> qPCR results </font></td>
      </tr>
      </table>
      <br/>
      <br/>
      <xsl:apply-templates select='PLATE'/>
    </body>

  </html>
</xsl:template>

<xsl:template match='PLATE'>
  <h2 align="center"> <xsl:value-of select="@NAME"/> </h2>
  <table cellpadding="2" cellspacing="0" border="1" width="100%" 
         bordercolor="black" borderstyle="solid">
      <tr>
        <th align="center" bgcolor="#ffffe0">Well</th>
        <th align="center" bgcolor="#ffffe0">Type</th>
        <th align="center" bgcolor="#ffffe0">Target</th>
        <th align="center" bgcolor="#ffffe0">Sample</th>
        <th align="center" bgcolor="#ffffe0">Ct</th>
        <th align="center" bgcolor="#ffffe0">Ctmean</th>
        <th align="center" bgcolor="#ffffe0">Ctdev</th>
        <th align="center" bgcolor="#ffffe0">Amount</th>
        <th align="center" bgcolor="#ffffe0">Efficiency</th>
        <th align="center" bgcolor="#ffffe0">NRQ</th>
        <th align="center" bgcolor="#ffffe0">NRQerror</th>
      </tr>
      <xsl:apply-templates select='WELL'/>
  </table>
  <br/>
</xsl:template>
 
<xsl:template match='WELL'>
  <tr>
    <th align="center"><font size="2"> <xsl:value-of select="NAME"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="TYPE"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="TARGET"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="SAMPLE"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="@CT"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="@CTMEAN"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="@CTDEV"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="@AMOUNT"/></font></th>
    <th align="center"><font size="2"> <xsl:apply-templates select='TARGET'/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="@NRQ"/></font></th>
    <th align="center"><font size="2"> <xsl:value-of select="@NRQERROR"/></font></th>
  </tr>
</xsl:template>

<xsl:template match='TARGET'>
  <xsl:value-of select="@EFF"/> &#177; <xsl:value-of select="@PM"/>
</xsl:template>

</xsl:stylesheet>
