import os
from pprint import pprint

# Update these paths according to your netMHCpan/netMHCIIpan downloads
NETMHCPAN40_PATH = '/Users/wsx/biosoft/netMHCpan-4.0/netMHCpan'
NETMHCIIPAN32_PATH = '/Users/wsx/biosoft/netMHCIIpan-3.2/netMHCIIpan'

# Update based on your desired tmp directory
TEMP_DIR = '/tmp/netMHC/'

if not os.path.exists(TEMP_DIR):
    os.makedirs(TEMP_DIR)

def setup_MHCI():
    return NETMHCPAN40_PATH, TEMP_DIR

def setup_MHCII():
    return NETMHCIIPAN32_PATH, TEMP_DIR


def check(obj):
    """
    Check all property for a object
    :param obj: a object
    :return: a pretty dictionary
    """
    pprint(vars(obj))