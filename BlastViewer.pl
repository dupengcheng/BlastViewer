#!/usr/bin/perl

=head1 Description
		draw linear schema of alignment between two sequences

=head1 Usage
		perl BlastViewer.pl [options]
			-ref      <file>             reference sequence file; necessary
			-qry      <file>             query sequence file; necessary
			-blast    <file>             blast parser file; necessary
			-ref_cds  <file>             reference cds sequence
			-qry_cds  <file>             query cds sequence
			-mark_cds <file>             list of marked cds
			-scale    <number>           scale for draw figure, default 2000
			-len_cut  <number>           blast length cutoff; default 1000
			-svg      <file>             output svg file; must be given
			-png                         create png file
			-help                        show help

=head1 Example

=head1 Version
		Author: Du Pengcheng, dupengcheng126@126.com
		Date: 2016-01-01

=cut


use Svg::File;
use Svg::Graphics;
use File::Basename;
use Getopt::Long;

my ($ref_file,$sca_file,$ref_cds_file,$qry_cds_file,$mark_cds,$blast_file,$scale,$filt,$svg_file);

GetOptions (
	"ref:s"=>\$ref_file,
	"qry:s"=>\$sca_file,
	"ref_cds:s"=>\$ref_cds_file,
	"qry_cds:s"=>\$qry_cds_file,
	"mark_cds:s"=>\$mark_cds,
	"blast:s"=>\$blast_file,
	"scale:i"=>\$scale,
	"len_cut:i"=>\$filt,
	"svg:s"=>\$svg_file,
	"png"=>\$png,
	"help"=>\$help,
);

die `pod2text $0` if ($help || !defined $blast_file || !defined $svg_file);

$scale=2000 unless ($scale);
$filt=1000 unless ($filt);

my $svg = Svg::File->new("$svg_file",\*OUT);
my $g = $svg->beginGraphics();
$svg->open(encoding, "iso-8859-1");
#$svg->open(silent);
$g->b(svg,viewBox,"0,0,2000,800",width,2000,height,800);


my $x0=200;
my $y0=200;
my $PI=3.1415;
my %color=("+"=>"#CC3366","-"=>"#3366FF");
my %cds_color=("+"=>"red","-"=>"green");

#############################################
#reference
#############################################
my (%ref,%sort_ref,%ref_sort,%ref_shift,%mark_cds,$ref_num);
if ($ref_file) {
	my $ref_name=basename $ref_file;
	$ref_name=~s/\.\w+//;
	$g->d(txtCM,"$ref_name",xval,$x0-100,yval,$y0,"doStyle");

	$/="\>";
	open (R,$ref_file) || die;
	<IN>;
	while (<R>) {
		my ($id, $seq)=split /\n/, $_, 2;
		$id=(split /\s/, $id)[0];
		$seq=~s/[\s\r\>\n]//g;
		next if (length($seq)==0);
		$ref{$id}=$seq;
	}


	close R;
	$/="\n";

	foreach (sort keys %ref) {
		$ref_num++;
		$sort_ref{$ref_num}=$_;
	}
	%ref_sort=reverse %sort_ref;

	foreach (keys %ref) {
		$ref_shift{$_}=0;
		for (my $i=1; $i<$ref_sort{$_}; $i++) {
			$ref_shift{$_}+=length($ref{$sort_ref{$i}})/$scale+10;
		}
	}

	foreach (keys %ref) {
		drawXRuler($x0+$ref_shift{$_},$y0,1,length($ref{$_}),$scale,100000,100000,"+","black");
	}
}

################################################################################
#scaffold
################################################################################
my (%sca,%sort_sca,%sca_sort,%sca_shift,$sca_num);
if ($sca_file) {
	my $sca_name=basename $sca_file;
	$sca_name=~s/\.\w+//;
	$g->d(txtCM,"$sca_name",xval,$x0-100,yval,$y0+200,"doStyle");

	my $id;
	open (SC,$sca_file) || die;
	while (<SC>) {
		chomp;
		if (/>(.[^\s]+)\s*/) {
			($id=$1)=~s/[\r\s]//g;
			$scaffold{$id}="";
			$sca_num++;
			$sort_sca{$sca_num}=$id;
		}
		else {
			s/\W//g;
			$sca{$id}.=uc($_);
		}
	}
	close SC;
	%sca_sort=reverse %sort_sca;

	foreach (keys %sca) {
		$sca_shift{$_}=0;
		for (my $i=1; $i<$sca_sort{$_}; $i++) {
			$sca_shift{$_}+=length($sca{$sort_sca{$i}})/$scale+10;
		}
	}

	foreach (keys %sca) {
		drawXRuler($x0+$sca_shift{$_},$y0+200,1,length($sca{$_}),$scale,100000,100000,"+","black");
	}
}

if ($mark_cds) {
	open (IN,$mark_cds) || die "OpenError: $mark_cds";
	while (<IN>) {
		chomp;
		my ($cds,$anno)=split;
		$mark_cds{$cds}=$anno;
	}
	close IN;
}

if ($ref_cds_file) {
	open (CDS1,"$ref_cds_file") || die "OpenError: $ref_cds_file";
	OUT:while (<CDS1>) {
		if (/>(\S+)\s+.+?[\:\=](.+?)\:(\d+)\:(\d+)\:([+-.])/) {
			my ($id, $seq, $str, $end, $dir)=($1,$2,$3,$4,$5);
			$seq=~s/[\s\r]//;
			next OUT if (!exists $ref{$seq});
			($str, $end)=sort {$a<=>$b} ($str, $end);
			my $color;
			if ($mark_cds) {
				$color=exists $mark_cds{$id} ? "blue":"gray";
				if (exists $mark_cds{$id}) {
					#$g->d(line,x1,($str+$end)/2/$scale+$x0+$ref_shift{$seq},y1,$y0-65,x2,($str+$end)/2/$scale+$x0+$ref_shift{$seq},y2,$y0-50,style,"stroke:$color;stroke-width:3");
					my $name_x=($str+$end)/2/$scale+$x0+$ref_shift{$seq};
					my $name_y=$y0-65;
					my $a=cos(-$PI/4);  my $b=sin(-$PI/4);  my $c=-sin(-$PI/4);  my $d=cos(-$PI/4);
					$g->d(txtLM,$mark_cds{$id},xval,0,yval,0,transform,"matrix($a $b $c $d $name_x $name_y)");
					#$g->d(txtLM, $mark_cds{$id}, xval, ($str+$end)/2/$scale+$x0+$ref_shift{$seq}, yval, $y0-65, "doStyle");
				}
			}
			else {
				$color=$cds_color{$dir};
			}
			drawXTriangleLine($str/$scale+$x0+$ref_shift{$seq},$end/$scale+$x0+$ref_shift{$seq},$y0-40,15,$dir,$color);
		}
	}
	close CDS1;
}

if ($qry_cds_file) {
	open (CDS2,"$qry_cds_file") || die "OpenError: $qry_cds_file";
	#print $qry_cds_file;
	OUT:while (<CDS2>) {
		if (/>(\S+)\s+.+?[\:\=](.+?)\:(\d+)\:(\d+)\:([+-.])/) {
			my ($id, $seq, $str, $end, $dir)=($1,$2,$3,$4,$5);
			next OUT if (!exists $sca{$seq});
			$seq=~s/[\s\r]//;
			($str, $end)=sort {$a<=>$b} ($str, $end);
			my $color;
			if ($mark_cds) {
				$color=exists $mark_cds{$id} ? "blue":"gray";
				if (exists $mark_cds{$id}) {
					my $name_x=($str+$end)/2/$scale+$x0+$sca_shift{$seq};
					my $name_y=$y0+235;
					my $a=cos($PI/4);  my $b=sin($PI/4);  my $c=-sin($PI/4);  my $d=cos($PI/4);
					$g->d(txtLM,$mark_cds{$id},xval,0,yval,0,transform,"matrix($a $b $c $d $name_x $name_y)");
					#$g->d(line,x1,($str+$end)/2/$scale+$x0+$sca_shift{$seq},y1,$y0+225,x2,($str+$end)/2/$scale+$x0+$sca_shift{$seq},y2,$y0+245,style,"stroke:$color;stroke-width:3");
					#$g->d(txtLM, $mark_cds{$id}, xval, ($str+$end)/2/$scale+$x0+$ref_shift{$seq}, yval, $y0+225, "doStyle");
				}
			}
			else {
				$color=$cds_color{$dir};
			}
			drawXTriangleLine($str/$scale+$x0+$sca_shift{$seq},$end/$scale+$x0+$sca_shift{$seq},$y0+220,15,$dir,$color);
		}
	}
	close CDS2;
}

open (BLA,"$blast_file") || die "OpenError: $blast_file";
while (<BLA>) {
	my ($sca,$sca_start,$sca_end,$ref,$ref_start,$ref_end,$idy,$align_len)=(split /\t/,$_)[0,2,3,4,6,7,8,11];
	if ($align_len<$filt) {
		next;
	}
	my $x1=$ref_start/$scale+$x0+$ref_shift{$ref}; my $y1=$y0+20;
	my $x2=$ref_end/$scale+$x0+$ref_shift{$ref}; my $y2=$y0+20;
	my $x3=$sca_end/$scale+$x0+$sca_shift{$sca}; my $y3=$y0+180;
	my $x4=$sca_start/$scale+$x0+$sca_shift{$sca}; my $y4=$y0+180;
	#my $color=$ref_start<$ref_end ? $color{"+"}:$color{"-"};
	my $gg=350*(1-$idy);
	#if ($idy<0.5) {
	#	$gg=200;
	#}
	#$g->d(polygon,style,"fill:$color;stroke:none;stoke-width:0",points,"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4");
	#$g->d(polygon,style,"fill:rgb(255,$gg,80);stroke:black;stoke-width:1",points,"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4");
	$g->d(polygon,style,"fill:rgb(0,$gg,255);stroke:none;stoke-width:0",points,"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4");
	$idy=$idy*100;
	$g->d(txtCM, $idy."%", xval, ($x1+$x2+$x3+$x4)/4, yval, ($y1+$y2+$y3+$y4)/4, "doStyle");
}


#$y0+=200;
#foreach (1..10) {
#	my $gg=350*(1-$_*0.1);
#	$g->d(rect,xval,$x0,yval,$y0+10*$_,width,20,height,10,style,"fill:rgb(0,$gg,255);stroke:none;stoke-width:0");
#}

#########################################################
#���������
$g->e();
$svg->close($g);	


#Usage: drawScale(xval,yval,scale_start,scale_end,scale,mini scale,label position,direction,color)
sub drawXRuler {
	my ($x,$y,$start,$end,$scale,$min_scale,$label_position,$direction,$color)=@_;

	#���������
	my $line_length=8;
	my $stroke_width=2;
	$g->d("line","x1",$x,"y1",$y,"x2",$x+($end-$start+1)/$scale,"y2",$y,style,"stroke:$color;stroke-width:$stroke_width");
	$g->d("line","x1",$x,"y1",$y,"x2",$x,"y2",$y-$line_length,style,"stroke:$color;stroke-width:$stroke_width");
	$g->d("line","x1",$x+($end-$start+1)/$scale,"y1",$y,"x2",$x+($end-$start+1)/$scale,"y2",$y-$line_length,style,"stroke:$color;stroke-width:$stroke_width");

	#�ж��Ƿ���Ҫ��ת
	my $reverse_flag;
	if ($direction eq "+") {
		$reverse_flag=1;
	}else {
		$x+=($end-$start+1)/$scale;
		$reverse_flag=-1;
	}

	for (my $i=$start; $i<$end; $i+=$min_scale) {
		if ($i==$start or $i%$label_position == 0) {
			$line_length=8;
			$stroke_width=2;
			$g->d("line","x1",$x+$reverse_flag*($i-$start)/$scale,"y1",$y-$line_length,"x2",$x+$reverse_flag*($i-$start)/$scale,"y2",$y,style,"stroke:$color;stroke-width:$stroke_width");

			if ($i==0) {
				$g->d(txtCM,$i+1,xval,$x+$reverse_flag*($i-$start)/$scale,yval,$y-$line_length-10);
			}
			elsif ($i==$start and $start%abs($min_scale)!=0) {
				$g->d(txtCM,$start,xval,$x+$reverse_flag*($i-$start)/$scale,yval,$y-$line_length-10);
				$i=int($start/abs($min_scale))*$min_scale;
			}
			else {
				my $value_k=($i/1000)."k";
				$g->d(txtCM,$value_k,xval,$x+$reverse_flag*($i-$start)/$scale,yval,$y-$line_length-10);
			}
		}
		else {
			$line_length=4;
			$stroke_width=1;
			$g->d("line","x1",$x+$reverse_flag*($i-$start)/$scale,"y1",$y-$line_length,"x2",$x+$reverse_flag*($i-$start)/$scale,"y2",$y,style,"stroke:$color;stroke-width:$stroke_width");
			#$g->e();
		}
	}

	if ($end%$label_position != 0) {
		$line_length=8;
		$stroke_width=2;
		$g->d("line","x1",$x+$reverse_flag*($end-$start)/$scale,"y1",$y-$line_length,"x2",$x+$reverse_flag*($end-$start)/$scale,"y2",$y,style,"stroke:$color;stroke-width:$stroke_width");
		#$g->e();
		if ($end%$min_scale != 0) {
			$g->d(txtCM,$end,xval,$x+$reverse_flag*($end-$start)/$scale,yval,$y-$line_length-10);
		} else {
			my $value_k=($end/1000)."k";
			$g->d(txtCM,$value_k,xval,$x+$reverse_flag*($end-$start)/$scale,yval,$y-$line_length-10);
		}
	}
}

####################################################################33
#Usage: drawArrowLine($x_start,$x_end,$y,$line_width,$direction,$color)
sub drawXTriangleLine {
	my ($x_start,$x_end,$y,$line_width,$direction,$color)=@_;
	my $arrow_length;
	$y-=0.5*$line_width;
	
	if ($direction eq ".") {
		$g->d(rect, xval, $x_start, yval, $y, width, $x_end-$x_start, height, $line_width, style,"fill:$color;stroke:none;stroke-width:0");
	}
	elsif ($x_end-$x_start > 0.75*1.732*$line_width and $direction eq "+") {
		$arrow_length=0.75*1.732*$line_width;
		#$g->d(rect,style,"fill:$color;stroke:black;stoke-width:0.5",xval,$x_start,yval,$y,width,$x_end-$x_start+1-$arrow_length,height,$line_width);
		my $rect_x1=$x_start;                  my $rect_y1=$y;
		my $rect_x2=$x_end-$arrow_length;      my $rect_y2=$y;
		my $rect_x3=$x_end-$arrow_length;      my $rect_y3=$y+$line_width;
		my $rect_x4=$x_start;                  my $rect_y4=$y+$line_width;
		my $triangle_x1=$x_end-$arrow_length;  my $triangle_y1=$y-$line_width/4;
		my $triangle_x2=$x_end;                my $triangle_y2=$y+$line_width/2;
		my $triangle_x3=$x_end-$arrow_length;  my $triangle_y3=$y+1.25*$line_width;
		$g->d(polygon,style,"fill:$color;stroke:none;stoke-width:0",points,"$rect_x1,$rect_y1 $rect_x2,$rect_y2 $triangle_x1,$triangle_y1 $triangle_x2,$triangle_y2 $triangle_x3,$triangle_y3 $rect_x3,$rect_y3 $rect_x4,$rect_y4");
		#$g->d(polygon,style,"fill:$color;stroke:black;stoke-width:0.5",points,"$triangle_x1,$triangle_y1 $triangle_x2,$triangle_y2 $triangle_x3,$triangle_y3");
	}
	elsif ($x_end-$x_start <= 0.75*1.732*$line_width and $direction eq "+") {
		$arrow_length=$x_end-$x_start+1;
		my $triangle_x1=$x_start;  my $triangle_y1=$y-$line_width/4;
		my $triangle_x2=$x_end;    my $triangle_y2=$y+$line_width/2;
		my $triangle_x3=$x_start;  my $triangle_y3=$y+1.25*$line_width;
		$g->d(polygon,style,"fill:$color;stroke:none;stoke-width:0",points,"$triangle_x1,$triangle_y1 $triangle_x2,$triangle_y2 $triangle_x3,$triangle_y3");
	}
	elsif ($x_end-$x_start > 0.75*1.732*$line_width and $direction eq "-") {
		$arrow_length=0.75*1.732*$line_width;
		#$g->d(rect,style,"fill:$color;stroke:black;stoke-width:1",xval,$x_start+$arrow_length,yval,$y,width,$x_end-$x_start+1-$arrow_length,height,$line_width);
		my $rect_x1=$x_end;                      my $rect_y1=$y;
		my $rect_x2=$x_start+$arrow_length;      my $rect_y2=$y;
		my $rect_x3=$x_start+$arrow_length;      my $rect_y3=$y+$line_width;
		my $rect_x4=$x_end;                      my $rect_y4=$y+$line_width;
		my $triangle_x1=$x_start+$arrow_length;  my $triangle_y1=$y-$line_width/4;
		my $triangle_x2=$x_start;                my $triangle_y2=$y+$line_width/2;
		my $triangle_x3=$x_start+$arrow_length;  my $triangle_y3=$y+1.25*$line_width;
		$g->d(polygon,style,"fill:$color;stroke:none;stoke-width:0",points,"$rect_x1,$rect_y1 $rect_x2,$rect_y2 $triangle_x1,$triangle_y1 $triangle_x2,$triangle_y2 $triangle_x3,$triangle_y3 $rect_x3,$rect_y3 $rect_x4,$rect_y4");
		#$g->d(polygon,style,"fill:$color;stroke:black;stoke-width:0.5",points,"$triangle_x1,$triangle_y1 $triangle_x2,$triangle_y2 $triangle_x3,$triangle_y3");
	}
	elsif ($x_end-$x_start <= 0.75*1.732*$line_width and $direction eq "-") {
		$arrow_length=$x_end-$x_start+1;
		my $triangle_x1=$x_end;      my $triangle_y1=$y-$line_width/4;
		my $triangle_x2=$x_start;    my $triangle_y2=$y+$line_width/2;
		my $triangle_x3=$x_end;      my $triangle_y3=$y+1.25*$line_width;
		$g->d(polygon,style,"fill:$color;stroke:none;stoke-width:0",points,"$triangle_x1,$triangle_y1 $triangle_x2,$triangle_y2 $triangle_x3,$triangle_y3");
	}
}
