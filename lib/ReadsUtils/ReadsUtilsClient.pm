package ReadsUtils::ReadsUtilsClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

ReadsUtils::ReadsUtilsClient

=head1 DESCRIPTION


Utilities for handling reads files.


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => ReadsUtils::ReadsUtilsClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = Bio::KBase::AuthToken->new(@args);
	
	if (!$token->error_message)
	{
	    $self->{token} = $token->token;
	    $self->{client}->{token} = $token->token;
	}
        else
        {
	    #
	    # All methods in this module require authentication. In this case, if we
	    # don't have a token, we can't continue.
	    #
	    die "Authentication failed: " . $token->error_message;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 validateFASTQ

  $out = $obj->validateFASTQ($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a list where each element is a ReadsUtils.ValidateFASTQParams
$out is a reference to a list where each element is a ReadsUtils.ValidateFASTQOutput
ValidateFASTQParams is a reference to a hash where the following keys are defined:
	file_path has a value which is a string
	interleaved has a value which is a ReadsUtils.boolean
boolean is an int
ValidateFASTQOutput is a reference to a hash where the following keys are defined:
	validated has a value which is a ReadsUtils.boolean

</pre>

=end html

=begin text

$params is a reference to a list where each element is a ReadsUtils.ValidateFASTQParams
$out is a reference to a list where each element is a ReadsUtils.ValidateFASTQOutput
ValidateFASTQParams is a reference to a hash where the following keys are defined:
	file_path has a value which is a string
	interleaved has a value which is a ReadsUtils.boolean
boolean is an int
ValidateFASTQOutput is a reference to a hash where the following keys are defined:
	validated has a value which is a ReadsUtils.boolean


=end text

=item Description

Validate a FASTQ file. The file extensions .fq, .fnq, and .fastq
are accepted. Note that prior to validation the file will be altered in
place to remove blank lines if any exist.

=back

=cut

 sub validateFASTQ
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function validateFASTQ (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to validateFASTQ:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'validateFASTQ');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "ReadsUtils.validateFASTQ",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'validateFASTQ',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method validateFASTQ",
					    status_line => $self->{client}->status_line,
					    method_name => 'validateFASTQ',
				       );
    }
}
 


=head2 upload_reads

  $return = $obj->upload_reads($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a ReadsUtils.UploadReadsParams
$return is a ReadsUtils.UploadReadsOutput
UploadReadsParams is a reference to a hash where the following keys are defined:
	fwd_id has a value which is a string
	wsid has a value which is an int
	wsname has a value which is a string
	objid has a value which is an int
	name has a value which is a string
	rev_id has a value which is a string
	sequencing_tech has a value which is a string
	single_genome has a value which is a ReadsUtils.boolean
	strain has a value which is a KBaseCommon.StrainInfo
	source has a value which is a KBaseCommon.SourceInfo
	interleaved has a value which is a ReadsUtils.boolean
	read_orientation_outward has a value which is a ReadsUtils.boolean
	insert_size_mean has a value which is a float
	insert_size_std_dev has a value which is a float
boolean is an int
StrainInfo is a reference to a hash where the following keys are defined:
	genetic_code has a value which is an int
	genus has a value which is a string
	species has a value which is a string
	strain has a value which is a string
	organelle has a value which is a string
	source has a value which is a KBaseCommon.SourceInfo
	ncbi_taxid has a value which is an int
	location has a value which is a KBaseCommon.Location
SourceInfo is a reference to a hash where the following keys are defined:
	source has a value which is a string
	source_id has a value which is a KBaseCommon.source_id
	project_id has a value which is a KBaseCommon.project_id
source_id is a string
project_id is a string
Location is a reference to a hash where the following keys are defined:
	lat has a value which is a float
	lon has a value which is a float
	elevation has a value which is a float
	date has a value which is a string
	description has a value which is a string
UploadReadsOutput is a reference to a hash where the following keys are defined:
	obj_ref has a value which is a string

</pre>

=end html

=begin text

$params is a ReadsUtils.UploadReadsParams
$return is a ReadsUtils.UploadReadsOutput
UploadReadsParams is a reference to a hash where the following keys are defined:
	fwd_id has a value which is a string
	wsid has a value which is an int
	wsname has a value which is a string
	objid has a value which is an int
	name has a value which is a string
	rev_id has a value which is a string
	sequencing_tech has a value which is a string
	single_genome has a value which is a ReadsUtils.boolean
	strain has a value which is a KBaseCommon.StrainInfo
	source has a value which is a KBaseCommon.SourceInfo
	interleaved has a value which is a ReadsUtils.boolean
	read_orientation_outward has a value which is a ReadsUtils.boolean
	insert_size_mean has a value which is a float
	insert_size_std_dev has a value which is a float
boolean is an int
StrainInfo is a reference to a hash where the following keys are defined:
	genetic_code has a value which is an int
	genus has a value which is a string
	species has a value which is a string
	strain has a value which is a string
	organelle has a value which is a string
	source has a value which is a KBaseCommon.SourceInfo
	ncbi_taxid has a value which is an int
	location has a value which is a KBaseCommon.Location
SourceInfo is a reference to a hash where the following keys are defined:
	source has a value which is a string
	source_id has a value which is a KBaseCommon.source_id
	project_id has a value which is a KBaseCommon.project_id
source_id is a string
project_id is a string
Location is a reference to a hash where the following keys are defined:
	lat has a value which is a float
	lon has a value which is a float
	elevation has a value which is a float
	date has a value which is a string
	description has a value which is a string
UploadReadsOutput is a reference to a hash where the following keys are defined:
	obj_ref has a value which is a string


=end text

=item Description

Loads a set of reads to KBase data stores.

=back

=cut

 sub upload_reads
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function upload_reads (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to upload_reads:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'upload_reads');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "ReadsUtils.upload_reads",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'upload_reads',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method upload_reads",
					    status_line => $self->{client}->status_line,
					    method_name => 'upload_reads',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "ReadsUtils.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "ReadsUtils.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'upload_reads',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method upload_reads",
            status_line => $self->{client}->status_line,
            method_name => 'upload_reads',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for ReadsUtils::ReadsUtilsClient\n";
    }
    if ($sMajor == 0) {
        warn "ReadsUtils::ReadsUtilsClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 boolean

=over 4



=item Description

A boolean - 0 for false, 1 for true.
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 ValidateFASTQParams

=over 4



=item Description

Input to the validateFASTQ function.

    Required parameters:
    file_path - the path to the file to validate.
    
    Optional parameters:
    interleaved - whether the file is interleaved or not. Setting this to
        true disables sequence ID checks.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
file_path has a value which is a string
interleaved has a value which is a ReadsUtils.boolean

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
file_path has a value which is a string
interleaved has a value which is a ReadsUtils.boolean


=end text

=back



=head2 ValidateFASTQOutput

=over 4



=item Description

The output of the validateFASTQ function.

validated - whether the file validated successfully or not.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
validated has a value which is a ReadsUtils.boolean

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
validated has a value which is a ReadsUtils.boolean


=end text

=back



=head2 UploadReadsParams

=over 4



=item Description

Input to the upload_reads function.

Required parameters:
fwd_id - the id of the shock node containing the reads data file:
    either single end reads, forward/left reads, or interleaved reads.
sequencing_tech - the sequencing technology used to produce the
    reads.

One of:
wsid - the id of the workspace where the reads will be saved
    (preferred).
wsname - the name of the workspace where the reads will be saved.

One of:
objid - the id of the workspace object to save over
name - the name to which the workspace object will be saved
    
Optional parameters:
rev_id - the shock node id containing the reverse/right reads for
    paired end, non-interleaved reads.
single_genome - whether the reads are from a single genome or a
    metagenome. Default is single genome.
strain - information about the organism strain
    that was sequenced.
source - information about the organism source.
interleaved - specify that the fwd reads file is an interleaved paired
    end reads file as opposed to a single end reads file. Default true,
    ignored if rev_id is specified.
read_orientation_outward - whether the read orientation is outward
    from the set of primers. Default is false and is ignored for
    single end reads.
insert_size_mean - the mean size of the genetic fragments. Ignored for
    single end reads.
insert_size_std_dev - the standard deviation of the size of the
    genetic fragments. Ignored for single end reads.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
fwd_id has a value which is a string
wsid has a value which is an int
wsname has a value which is a string
objid has a value which is an int
name has a value which is a string
rev_id has a value which is a string
sequencing_tech has a value which is a string
single_genome has a value which is a ReadsUtils.boolean
strain has a value which is a KBaseCommon.StrainInfo
source has a value which is a KBaseCommon.SourceInfo
interleaved has a value which is a ReadsUtils.boolean
read_orientation_outward has a value which is a ReadsUtils.boolean
insert_size_mean has a value which is a float
insert_size_std_dev has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
fwd_id has a value which is a string
wsid has a value which is an int
wsname has a value which is a string
objid has a value which is an int
name has a value which is a string
rev_id has a value which is a string
sequencing_tech has a value which is a string
single_genome has a value which is a ReadsUtils.boolean
strain has a value which is a KBaseCommon.StrainInfo
source has a value which is a KBaseCommon.SourceInfo
interleaved has a value which is a ReadsUtils.boolean
read_orientation_outward has a value which is a ReadsUtils.boolean
insert_size_mean has a value which is a float
insert_size_std_dev has a value which is a float


=end text

=back



=head2 UploadReadsOutput

=over 4



=item Description

The output of the upload_reads function.

    obj_ref - a reference to the new Workspace object in the form X/Y/Z,
        where X is the workspace ID, Y is the object ID, and Z is the
        version.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
obj_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
obj_ref has a value which is a string


=end text

=back



=cut

package ReadsUtils::ReadsUtilsClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
