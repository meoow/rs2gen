package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"gomiscutils"
	"log"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
)

const (
	NULL uint8 = iota
	P3UPSTREAM
	P5UPSTREAM
	P3UTR
	P5UTR
	INTRON
	CODING
)
const VERSION string = "0.2.0a"

type NullValue struct{}

type geneInfo struct {
	geneType uint8
	geneName map[string]*geneInfo_
}
type geneInfo_ struct {
	geneID   uint32
	geneCode uint8
}

type rsID uint32
type rsMap map[rsID]NullValue
type rsInfo map[rsID]*geneInfo

var dbsnpfile string
var snplist string
var norule bool
var version bool

func init() {
	flag.StringVar(&dbsnpfile, "db", "DBSNP.GZ", "Database file from dbSNP (gz compressed)")
	flag.StringVar(&snplist, "rs", "RSLIST.TXT", "One rs number per line listed file.")
	flag.BoolVar(&norule,"nr",false,"Do not map according to the gene priority rules.")
	flag.BoolVar(&version, "v", false, "Show version")
}

func main() {
	flag.Parse()
	if version {
		fmt.Println(VERSION)
		os.Exit(1)
	}
//	if flag.NFlag() != 2 {
//		flag.PrintDefaults()
//		os.Exit(0)
//	}
	b138file := dbsnpfile
	includeLiseFile := snplist
	includeList := readrslist(includeLiseFile)
	b138List := parseB138data(b138file, includeList)
	rs2gene(b138List)
}

func die(e error) {
	logger := log.New(os.Stderr, "[SHIT!!] ", 0)
	if e != nil {
		logger.Fatal(e)
	}
}

func readrslist(pathOfFile string) rsMap {
	var filer *os.File
	var err error
	var fReader *bufio.Reader
	var rsList rsMap = make(rsMap)

	filer, err = os.Open(pathOfFile)
	die(err)
	defer filer.Close()
	fReader = bufio.NewReader(filer)
	lines := gomiscutils.Readline(fReader)
	for l := range lines {
		if strings.HasPrefix(l, "rs") {
			l = l[2:]
		}
		k, err := strconv.ParseUint(gomiscutils.TrimNewLine(l), 10, 32)
		die(err)
		rsList[rsID(k)] = NullValue{}
	}
	return rsList
}

func parseB138data(pathOfFile string, includeList rsMap) rsInfo {
	var filer *os.File
	var err error
	var fReader *bufio.Reader
	var gzipReader *gzip.Reader
	var b138RS rsID
	var b138RSInfo = make(rsInfo, len(includeList))
	var linecount uint32
	var linecountStr string
	var linecountStrLen int

	rand.Seed(time.Now().Unix())
	var randLines = rand.Uint32() / 100000

	filer, err = os.Open(pathOfFile)
	die(err)
	defer filer.Close()

	gzipReader, err = gzip.NewReader(filer)
	die(err)

	fReader = bufio.NewReader(gzipReader)
	lines := gomiscutils.Readline(fReader)
	os.Stderr.WriteString("Scanning SNP database...0")
	linecountStrLen = 1
	for l := range lines {
		linecount++
		if linecount%randLines == 0 {
			randLines = rand.Uint32() % 100000
			os.Stderr.WriteString(strings.Repeat("\b", linecountStrLen))
			linecountStr = fmt.Sprintf("%d", linecount)
			linecountStrLen = len(linecountStr)
			os.Stderr.WriteString(linecountStr)
		}
		//fmt.Printf("%v",l)
		//l = strings.TrimRight(l, newLineChar)
		fields := strings.Split(l, "\t")
		if fields[5] == "" || fields[6] == "" || fields[1][0:3] == "NW_" {
			continue
		}
		//b138RS map[rsID]map[geneName]*geneInfo
		b138RS_, err := strconv.ParseUint(fields[0], 10, 32)
		die(err)
		b138RS = rsID(b138RS_)
		if _, ok := includeList[b138RS]; ok {
			geneid_, err := strconv.ParseUint(fields[5], 10, 32)
			die(err)
			geneid := uint32(geneid_)
			genename := fields[6]
			genecode_, err := strconv.ParseUint(fields[11], 10, 8)
			genecode := uint8(genecode_)
			if _, ok := b138RSInfo[b138RS]; !ok {
				b138RSInfo[b138RS] = &geneInfo{NULL, make(map[string]*geneInfo_)}
			}
			b138RSInfo[b138RS].addGene(genecode, genename, geneid)
		}

	}
	os.Stderr.WriteString(strings.Repeat("\b", linecountStrLen))
	linecountStr = fmt.Sprintf("%d", linecount)
	os.Stderr.WriteString(linecountStr)
	os.Stderr.WriteString("\n")
	return b138RSInfo
}

func rs2gene(r rsInfo) {
	for rs, gi := range r {
		for gn, gi_ := range gi.geneName {
			fmt.Printf("%d\t%d\t%s\t%d\n", rs, gi_.geneID, gn, gi_.geneCode)
		}
	}
}

func (gi *geneInfo) addGene(genecode uint8, genename string, geneid uint32) {
	if !norule {
		var setInfo func(code uint8) = func(code uint8) {
			if code > gi.geneType {
				gi.geneType = code
				for k, _ := range gi.geneName {
					delete(gi.geneName, k)
				}
				gi.geneName[genename] = &geneInfo_{geneid, genecode}
			} else if code == gi.geneType {
				if _, ok := gi.geneName[genename]; ok {
					return
				} else {
					gi.geneName[genename] = &geneInfo_{geneid, genecode}
				}
			}
		}
		switch genecode {
		case 20, 3, 8, 9, 41, 42, 43, 44, 45, 30:
			setInfo(CODING)
		case 6, 75, 73:
			setInfo(INTRON)
		case 55:
			setInfo(P5UTR)
		case 53:
			setInfo(P3UTR)
		case 15:
			setInfo(P5UPSTREAM)
		case 13:
			setInfo(P3UPSTREAM)
		}
		return
	} else {
		gi.geneName[genename] = &geneInfo_{geneid, genecode}
	}
}
